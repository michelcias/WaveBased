#include <R.h>
#include <Rinternals.h>

/* This function calculates the minimum (xmin) and the maximum (xmax)
 * among all the finite observations of the data x
 */
void Range(double *xmin, double *xmax, double *x, int n){
   int i, k = 0;
   // Initializing xmin as the first finite observation
   for(i = 0; i < n; i++){
     if(R_FINITE(x[i])){
       xmin[0] = x[i];
       k = i;
       break;
     }
   }
   
   // (Just being careful!) If there is no finite observation
   // in the data x, finish the code
   if(ISNA(xmin[0]))
     return;
   
   // Initializing xmax as the first finite observation
   xmax[0] = xmin[0];
   
   // Selecting the minimum (xmin) and the maximum (xmax)
   // among all the finite observations
   for(i = k; i < n; i++){
     if(xmin[0] > x[i] && R_FINITE(x[i]))
       xmin[0] = x[i];
     
     if(xmax[0] < x[i] && R_FINITE(x[i]))
       xmax[0] = x[i];
   }
   
 }

/* This function is useful to periodize the auxiliary matrix 'amat' and set the
 * periodized matrix as 'pmat'
 */
void Periodize(double *pmat, double *amat, int nrows, int k1, int k2, int p){
  
  int i, j, k, l;
  
  for(i = 0; i < nrows; i++){
    if(!R_FINITE(amat[i])){
      for(k = 0; k < p; k++)
        pmat[i + nrows*k] = NA_REAL;
      continue;
    }
    
    for(j = k1; j < (k1 + p); j++){
      k = (j % p + p) % p;
      pmat[i + nrows*k] = 0.0;
      for(l = j; l < (k2 + 1); l += p)
        pmat[i + nrows*k] += amat[i + nrows*(l - k1)];
    }
  }
  
}

/* First of all, this function is based on the 'phi' fucntion in the file 
 * 'WAVDE.c' of the wavethresh packege. Basically, this 'PhiVec' calculates
 * phi[0k](x) for all k value in which phi is non-null.
 */
void PhiVec(double *phi, double px, double *filter, int N, int prec, double *prod, double *tmp){
  double tval, w;
  int d, i, ind, j, k, l;
  
  // ----- Calculating phi[0k](px) ----- //
  // Initializing prod as an identity matrix
  for(i = 0; i < (N - 1); i++)
    for(j = 0; j < (N - 1); j++)
      prod[i + (N - 1)*j] = (i == j) ? 1.0 : 0.0;
  
  w = px - floor(px);
  
  if(w == 0)
    w = w + 1e-9;
  else if(w == 1)
    w = w - 1e-9;
  
  for(l = 0; l < prec; l++){
    
    // Calculating the l-th coefficient of the diadic expansion of w
    // {d[l] --> d}
    w = 2*w;
    d = floor(w);
    
    for(i = 0; i < (N - 1); i++){
      for(j = 0; j < (N - 1); j++){
        tmp[i + (N - 1)*j] = 0;
        
        for(k = 0; k < (N - 1); k++){
          
          // Calculating {T[k,j] --> tval}
          ind = 2*(k + 1) - (j + 1) - 1 + d;
          if(ind < 0 || ind > (N - 1))
            tval = 0;
          else
            tval = sqrt(2)*filter[ind];
          
          tmp[i + (N - 1)*j] += prod[i + (N - 1)*k]*tval;
        }
      }
    }
    
    for(i = 0; i < (N - 1)*(N - 1); i++)
      prod[i] = tmp[i];
    
    w = w - d; // updating w to the next diadic coefficient
  }
  
  // Calculating phi[i] as the average of the columns from the prod matrix
  for(i = 0; i < (N - 1); i++){
    phi[N - 2 - i] = 0;
    for(j = 0; j < (N - 1); j++){
      phi[N - 2 - i] += prod[i + (N - 1)*j]/(N - 1);
    }
  }
  // ----- phi[0k](px) calculated  ----- //
}

/* 'PsiVec' calculates psi[0k](x) for all k value in which psi is non-null.
 * Here we consider the result in Theorem 3.5.5 of Vidakovic (2002, p. 90).
 */
void PsiVec(double *psi, double px, double *filter, int N, int prec, int kmin, double *prod, double *tmp){
  
  double px2 = 2*px, uval, *phi;
  int i, j, fpx2;
  int *idxu;
  
  fpx2 = floor(px2);
  idxu = (int *) R_alloc((3*(N - 2) + 1), sizeof(int));
  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  
  PhiVec(phi, px2, filter, N, prec, prod, tmp);
  
  for(i = 0; i < (3*(N - 2) + 1); i++)
    idxu[i] = i + 1 + 2 * kmin - fpx2;
  
  // ----- Calculating psi[0k](px) based on {psi <-- Umat %*% phi} ----- //
  for(i = 0; i < (N - 1); i++){
    psi[i] = 0;
    for(j = 0; j < (N - 1); j++){
      if(idxu[2*i + j] < 0 || idxu[2*i + j] > (N - 1))
        uval = 0.0;
      else
        uval = pow((-1), idxu[2*i + j] - 1) * filter[idxu[2*i + j]];
      
      psi[i] += sqrt(2) * uval * phi[(N - 1) - j - 1];
    }
  }
  // ----- psi[0k](px) calculated  ----- //
}

/* For the sake of computational efficiency, this function
 * provides the filter and the limit points of the support of phi and psi, based
 * on the family ("daublets" [1], "symmlets" [2] or "coiflets" [3]) and the
 * number ot taps of the filter. All these filters were taken from the
 * function 'filter.select' of the wavethresh package.
 */
SEXP WavUtilities(SEXP family, SEXP taps){
  
  double *rwfilter;
  int rfam, rtaps;
  rfam = INTEGER(family)[0];
  rtaps = INTEGER(taps)[0];
  
  SEXP wfilter = PROTECT(allocVector(REALSXP, rtaps));
  rwfilter = REAL(wfilter);
  
  if(rfam == 1){
    if(rtaps == 2){
      rwfilter[0] = 1/sqrt(2);
      rwfilter[1] = 1/sqrt(2);
    }
    else if(rtaps == 4){
      rwfilter[0] =  0.48296291314469025;
      rwfilter[1] =  0.836516303737469;
      rwfilter[2] =  0.22414386804185735;
      rwfilter[3] = -0.12940952255092145;
    }
    else if(rtaps == 6){
      rwfilter[0] =  0.3326705529509569;
      rwfilter[1] =  0.8068915093133388;
      rwfilter[2] =  0.4598775021193313;
      rwfilter[3] = -0.13501102001039084;
      rwfilter[4] = -0.08544127388224149;
      rwfilter[5] =  0.035226291882100656;
    }
    else if(rtaps == 8){
      rwfilter[0] =  0.23037781330885523;
      rwfilter[1] =  0.7148465705525415;
      rwfilter[2] =  0.6308807679295904;
      rwfilter[3] = -0.02798376941698385;
      rwfilter[4] = -0.18703481171888114;
      rwfilter[5] =  0.030841381835986965;
      rwfilter[6] =  0.032883011666982945;
      rwfilter[7] = -0.010597401784997278;
    }
    else if(rtaps == 10){
      rwfilter[0] =  0.160102397974125;
      rwfilter[1] =  0.6038292697974729;
      rwfilter[2] =  0.7243085284385744;
      rwfilter[3] =  0.13842814590110342;
      rwfilter[4] = -0.24229488706619015;
      rwfilter[5] = -0.03224486958502952;
      rwfilter[6] =  0.07757149384006515;
      rwfilter[7] = -0.006241490213011705;
      rwfilter[8] = -0.012580751999015526;
      rwfilter[9] =  0.003335725285001549;
    }
    else if(rtaps == 12){
      rwfilter[0]  =  0.11154074335;
      rwfilter[1]  =  0.494623890398;
      rwfilter[2]  =  0.751133908021;
      rwfilter[3]  =  0.315250351709;
      rwfilter[4]  = -0.226264693965;
      rwfilter[5]  = -0.129766867567;
      rwfilter[6]  =  0.097501605587;
      rwfilter[7]  =  0.02752286553;
      rwfilter[8]  = -0.031582039318;
      rwfilter[9]  =  0.000553842201;
      rwfilter[10] =  0.004777257511;
      rwfilter[11] = -0.001077301085;
    }
    else if(rtaps == 14){
      rwfilter[0]  =  0.077852054085;
      rwfilter[1]  =  0.396539319482;
      rwfilter[2]  =  0.729132090846;
      rwfilter[3]  =  0.469782287405;
      rwfilter[4]  = -0.143906003929;
      rwfilter[5]  = -0.224036184994;
      rwfilter[6]  =  0.071309219267;
      rwfilter[7]  =  0.080612609151;
      rwfilter[8]  = -0.038029936935;
      rwfilter[9]  = -0.016574541631;
      rwfilter[10] =  0.012550998556;
      rwfilter[11] =  0.000429577973;
      rwfilter[12] = -0.001801640704;
      rwfilter[13] =  0.0003537138;
    }
    else if(rtaps == 16){
      rwfilter[0]  =  0.054415842243;
      rwfilter[1]  =  0.312871590914;
      rwfilter[2]  =  0.675630736297;
      rwfilter[3]  =  0.585354683654;
      rwfilter[4]  = -0.015829105256;
      rwfilter[5]  = -0.284015542962;
      rwfilter[6]  =  0.000472484574;
      rwfilter[7]  =  0.12874742662;
      rwfilter[8]  = -0.017369301002;
      rwfilter[9]  = -0.044088253931;
      rwfilter[10] =  0.013981027917;
      rwfilter[11] =  0.008746094047;
      rwfilter[12] = -0.004870352993;
      rwfilter[13] = -0.000391740373;
      rwfilter[14] =  0.000675449406;
      rwfilter[15] = -0.000117476784;
    }
    else if(rtaps == 18){
      rwfilter[0]  =  0.038077947364;
      rwfilter[1]  =  0.243834674613;
      rwfilter[2]  =  0.60482312369;
      rwfilter[3]  =  0.657288078051;
      rwfilter[4]  =  0.133197385825;
      rwfilter[5]  = -0.293273783279;
      rwfilter[6]  = -0.096840783223;
      rwfilter[7]  =  0.148540749338;
      rwfilter[8]  =  0.030725681479;
      rwfilter[9]  = -0.067632829061;
      rwfilter[10] =  0.000250947115;
      rwfilter[11] =  0.022361662124;
      rwfilter[12] = -0.004723204758;
      rwfilter[13] = -0.004281503682;
      rwfilter[14] =  0.001847646883;
      rwfilter[15] =  0.000230385764;
      rwfilter[16] = -0.000251963189;
      rwfilter[17] =  3.934732e-05;
    }
    else if(rtaps == 20){
      rwfilter[0]  =  0.026670057901;
      rwfilter[1]  =  0.188176800078;
      rwfilter[2]  =  0.527201188932;
      rwfilter[3]  =  0.688459039454;
      rwfilter[4]  =  0.281172343661;
      rwfilter[5]  = -0.249846424327;
      rwfilter[6]  = -0.195946274377;
      rwfilter[7]  =  0.127369340336;
      rwfilter[8]  =  0.093057364604;
      rwfilter[9]  = -0.071394147166;
      rwfilter[10] = -0.029457536822;
      rwfilter[11] =  0.033212674059;
      rwfilter[12] =  0.003606553567;
      rwfilter[13] = -0.010733175483;
      rwfilter[14] =  0.001395351747;
      rwfilter[15] =  0.001992405295;
      rwfilter[16] = -0.000685856695;
      rwfilter[17] = -0.000116466855;
      rwfilter[18] =  9.358867e-05;
      rwfilter[19] = -1.3264203e-05;
    }
    else
      error("'taps = %d' is not allowed for 'Daublets'. For this family, only 2, 4, 6, ..., 18 and 20 taps are available.", rtaps);
  }
  else if(rfam == 2){
    if(rtaps == 8){
      rwfilter[0] = -0.107148901418/sqrt(2);
      rwfilter[1] = -0.041910965125/sqrt(2);
      rwfilter[2] =  0.703739068656/sqrt(2);
      rwfilter[3] =  1.136658243408/sqrt(2);
      rwfilter[4] =  0.421234534204/sqrt(2);
      rwfilter[5] = -0.140317624179/sqrt(2);
      rwfilter[6] = -0.017824701442/sqrt(2);
      rwfilter[7] =  0.045570345896/sqrt(2);
    }
    else if(rtaps == 10){
      rwfilter[0]  =  0.038654795955/sqrt(2);
      rwfilter[1]  =  0.041746864422/sqrt(2);
      rwfilter[2]  = -0.055344186117/sqrt(2);
      rwfilter[3]  =  0.281990696854/sqrt(2);
      rwfilter[4]  =  1.023052966894/sqrt(2);
      rwfilter[5]  =  0.89658164838/sqrt(2);
      rwfilter[6]  =  0.023478923136/sqrt(2);
      rwfilter[7]  = -0.247951362613/sqrt(2);
      rwfilter[8]  = -0.029842499869/sqrt(2);
      rwfilter[9]  =  0.027632152958/sqrt(2);
    }
    else if(rtaps == 12){
      rwfilter[0]  =  0.021784700327/sqrt(2);
      rwfilter[1]  =  0.004936612372/sqrt(2);
      rwfilter[2]  = -0.166863215412/sqrt(2);
      rwfilter[3]  = -0.068323121587/sqrt(2);
      rwfilter[4]  =  0.694457972958/sqrt(2);
      rwfilter[5]  =  1.113892783926/sqrt(2);
      rwfilter[6]  =  0.477904371333/sqrt(2);
      rwfilter[7]  = -0.102724969862/sqrt(2);
      rwfilter[8]  = -0.029783751299/sqrt(2);
      rwfilter[9]  =  0.06325056266/sqrt(2);
      rwfilter[10] =  0.002499922093/sqrt(2);
      rwfilter[11] = -0.011031867509/sqrt(2);
    }
    else if(rtaps == 14){
      rwfilter[0]  =  0.003792658534/sqrt(2);
      rwfilter[1]  = -0.001481225915/sqrt(2);
      rwfilter[2]  = -0.017870431651/sqrt(2);
      rwfilter[3]  =  0.043155452582/sqrt(2);
      rwfilter[4]  =  0.096014767936/sqrt(2);
      rwfilter[5]  = -0.070078291222/sqrt(2);
      rwfilter[6]  =  0.024665659489/sqrt(2);
      rwfilter[7]  =  0.758162601964/sqrt(2);
      rwfilter[8]  =  1.085782709814/sqrt(2);
      rwfilter[9]  =  0.408183939725/sqrt(2);
      rwfilter[10] = -0.198056706807/sqrt(2);
      rwfilter[11] = -0.152463871896/sqrt(2);
      rwfilter[12] =  0.005671342686/sqrt(2);
      rwfilter[13] =  0.014521394762/sqrt(2);
    }
    else if(rtaps == 16){
      rwfilter[0]  =  0.002672793393/sqrt(2);
      rwfilter[1]  = -0.0004283943/sqrt(2);
      rwfilter[2]  = -0.021145686528/sqrt(2);
      rwfilter[3]  =  0.005386388754/sqrt(2);
      rwfilter[4]  =  0.069490465911/sqrt(2);
      rwfilter[5]  = -0.038493521263/sqrt(2);
      rwfilter[6]  = -0.073462508761/sqrt(2);
      rwfilter[7]  =  0.515398670374/sqrt(2);
      rwfilter[8]  =  1.099106630537/sqrt(2);
      rwfilter[9]  =  0.68074534719/sqrt(2);
      rwfilter[10] = -0.086653615406/sqrt(2);
      rwfilter[11] = -0.202648655286/sqrt(2);
      rwfilter[12] =  0.010758611751/sqrt(2);
      rwfilter[13] =  0.044823623042/sqrt(2);
      rwfilter[14] = -0.000766690896/sqrt(2);
      rwfilter[15] = -0.004783458512/sqrt(2);
    }
    else if(rtaps == 18){
      rwfilter[0]  =  0.001512487309/sqrt(2);
      rwfilter[1]  = -0.000669141509/sqrt(2);
      rwfilter[2]  = -0.014515578553/sqrt(2);
      rwfilter[3]  =  0.012528896242/sqrt(2);
      rwfilter[4]  =  0.087791251554/sqrt(2);
      rwfilter[5]  = -0.02578644593/sqrt(2);
      rwfilter[6]  = -0.270893783503/sqrt(2);
      rwfilter[7]  =  0.049882830959/sqrt(2);
      rwfilter[8]  =  0.873048407349/sqrt(2);
      rwfilter[9]  =  1.015259790832/sqrt(2);
      rwfilter[10] =  0.337658923602/sqrt(2);
      rwfilter[11] = -0.077172161097/sqrt(2);
      rwfilter[12] =  0.000825140929/sqrt(2);
      rwfilter[13] =  0.042744433602/sqrt(2);
      rwfilter[14] = -0.016303351226/sqrt(2);
      rwfilter[15] = -0.018769396836/sqrt(2);
      rwfilter[16] =  0.000876502539/sqrt(2);
      rwfilter[17] =  0.001981193736/sqrt(2);
    }
    else if(rtaps == 20){
      rwfilter[0]  =  0.001089170447/sqrt(2);
      rwfilter[1]  =  0.00013524502/sqrt(2);
      rwfilter[2]  = -0.01222064263/sqrt(2);
      rwfilter[3]  = -0.002072363923/sqrt(2);
      rwfilter[4]  =  0.064950924579/sqrt(2);
      rwfilter[5]  =  0.016418869426/sqrt(2);
      rwfilter[6]  = -0.225558972234/sqrt(2);
      rwfilter[7]  = -0.100240215031/sqrt(2);
      rwfilter[8]  =  0.667071338154/sqrt(2);
      rwfilter[9]  =  1.0882515305/sqrt(2);
      rwfilter[10] =  0.542813011213/sqrt(2);
      rwfilter[11] = -0.050256540092/sqrt(2);
      rwfilter[12] = -0.045240772218/sqrt(2);
      rwfilter[13] =  0.07070356755/sqrt(2);
      rwfilter[14] =  0.008152816799/sqrt(2);
      rwfilter[15] = -0.028786231926/sqrt(2);
      rwfilter[16] = -0.001137535314/sqrt(2);
      rwfilter[17] =  0.006495728375/sqrt(2);
      rwfilter[18] =  8.0661204e-05/sqrt(2);
      rwfilter[19] = -0.000649589896/sqrt(2);
    }
    else
      error("'taps = %d' is not allowed for 'Symmlets'. For this family, only 10, 12, 14, 16, 18 and 20 taps are available.", rtaps);
  }
  else if(rfam == 3){
    /* I know. Some (many!) indices below should be negative here, but this
    * will not affect the final matrix basis of the Coiflets coefficients.
    * At the end of the day, the translation values of the current wavelet
    * basis will be shifted, and everything will be alright! ;-)
    */ 
    if(rtaps == 6){
      rwfilter[0] = -0.051429728471*sqrt(2);
      rwfilter[1] =  0.238929728471*sqrt(2);
      rwfilter[2] =  0.602859456942*sqrt(2);
      rwfilter[3] =  0.272140543058*sqrt(2);
      rwfilter[4] = -0.051429972847*sqrt(2);
      rwfilter[5] = -0.011070271529*sqrt(2);
    }
    else if(rtaps == 12){
      rwfilter[0]  =  0.0115876*sqrt(2);
      rwfilter[1]  = -0.02932014*sqrt(2);
      rwfilter[2]  = -0.04763959*sqrt(2);
      rwfilter[3]  =  0.273021*sqrt(2);
      rwfilter[4]  =  0.5746824*sqrt(2);
      rwfilter[5]  =  0.2948672*sqrt(2);
      rwfilter[6]  = -0.05408561*sqrt(2);
      rwfilter[7]  = -0.04202648*sqrt(2);
      rwfilter[8]  =  0.01674441*sqrt(2);
      rwfilter[9]  =  0.003967884*sqrt(2);
      rwfilter[10] = -0.001289203*sqrt(2);
      rwfilter[11] = -0.0005095054*sqrt(2);
    }
    else if(rtaps == 18){
      rwfilter[0]  = -0.002682419*sqrt(2);
      rwfilter[1]  =  0.005503127*sqrt(2);
      rwfilter[2]  =  0.01658356*sqrt(2);
      rwfilter[3]  = -0.04650776*sqrt(2);
      rwfilter[4]  = -0.04322076*sqrt(2);
      rwfilter[5]  =  0.2865033*sqrt(2);
      rwfilter[6]  =  0.5612853*sqrt(2);
      rwfilter[7]  =  0.3029836*sqrt(2);
      rwfilter[8]  = -0.05077014*sqrt(2);
      rwfilter[9]  = -0.05819625*sqrt(2);
      rwfilter[10] =  0.02443409*sqrt(2);
      rwfilter[11] =  0.01122924*sqrt(2);
      rwfilter[12] = -0.006369601*sqrt(2);
      rwfilter[13] = -0.001820459*sqrt(2);
      rwfilter[14] =  0.0007902051*sqrt(2);
      rwfilter[15] =  0.0003296652*sqrt(2);
      rwfilter[16] = -5.019277e-05*sqrt(2);
      rwfilter[17] = -2.446573e-05*sqrt(2);
    }
    else if(rtaps == 24){
      rwfilter[0]  =  0.000630961*sqrt(2);
      rwfilter[1]  = -0.001152225*sqrt(2);
      rwfilter[2]  = -0.005194524*sqrt(2);
      rwfilter[3]  =  0.01136246*sqrt(2);
      rwfilter[4]  =  0.01886724*sqrt(2);
      rwfilter[5]  = -0.05746423*sqrt(2);
      rwfilter[6]  = -0.03965265*sqrt(2);
      rwfilter[7]  =  0.2936674*sqrt(2);
      rwfilter[8]  =  0.5531265*sqrt(2);
      rwfilter[9]  =  0.3071573*sqrt(2);
      rwfilter[10] = -0.04711274*sqrt(2);
      rwfilter[11] = -0.06803813*sqrt(2);
      rwfilter[12] =  0.02781364*sqrt(2);
      rwfilter[13] =  0.01773584*sqrt(2);
      rwfilter[14] = -0.01075632*sqrt(2);
      rwfilter[15] = -0.004001013*sqrt(2);
      rwfilter[16] =  0.002652666*sqrt(2);
      rwfilter[17] =  0.0008955945*sqrt(2);
      rwfilter[18] = -0.0004165006*sqrt(2);
      rwfilter[19] = -0.0001838298*sqrt(2);
      rwfilter[20] =  4.408035e-05*sqrt(2);
      rwfilter[21] =  2.208286e-05*sqrt(2);
      rwfilter[22] = -2.304942e-06*sqrt(2);
      rwfilter[23] = -1.262175e-06*sqrt(2);
    }
    else if(rtaps == 30){
      rwfilter[0]  = -0.0001499638*sqrt(2);
      rwfilter[1]  =  0.0002535612*sqrt(2);
      rwfilter[2]  =  0.001540246*sqrt(2);
      rwfilter[3]  = -0.002941111*sqrt(2);
      rwfilter[4]  = -0.007163782*sqrt(2);
      rwfilter[5]  =  0.01655207*sqrt(2);
      rwfilter[6]  =  0.0199178*sqrt(2);
      rwfilter[7]  = -0.06499726*sqrt(2);
      rwfilter[8]  = -0.03680007*sqrt(2);
      rwfilter[9]  =  0.2980923*sqrt(2);
      rwfilter[10] =  0.5475054*sqrt(2);
      rwfilter[11] =  0.3097068*sqrt(2);
      rwfilter[12] = -0.04386605*sqrt(2);
      rwfilter[13] = -0.07465224*sqrt(2);
      rwfilter[14] =  0.02919588*sqrt(2);
      rwfilter[15] =  0.02311078*sqrt(2);
      rwfilter[16] = -0.01397369*sqrt(2);
      rwfilter[17] = -0.00648009*sqrt(2);
      rwfilter[18] =  0.004783001*sqrt(2);
      rwfilter[19] =  0.001720655*sqrt(2);
      rwfilter[20] = -0.001175822*sqrt(2);
      rwfilter[21] = -0.000451227*sqrt(2);
      rwfilter[22] =  0.0002137298*sqrt(2);
      rwfilter[23] =  9.93776e-05*sqrt(2);
      rwfilter[24] = -2.92321e-05*sqrt(2);
      rwfilter[25] = -1.5072e-05*sqrt(2);
      rwfilter[26] =  2.6408e-06*sqrt(2);
      rwfilter[27] =  1.4593e-06*sqrt(2);
      rwfilter[28] = -1.184e-07*sqrt(2);
      rwfilter[29] = -6.73e-08*sqrt(2);
    }
    else
      error("'taps = %d' is not allowed for 'Coiflets'. For this family, only 12, 18, 24 and 30 taps are available.", rtaps);
  }
  else
    error("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'.", rfam);
  
  SEXP results = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(results, 0, ScalarInteger(0));
  SET_VECTOR_ELT(results, 1, ScalarInteger(rtaps - 1));
  SET_VECTOR_ELT(results, 2, ScalarInteger(-(rtaps/2 - 1)));
  SET_VECTOR_ELT(results, 3, ScalarInteger(rtaps/2));
  SET_VECTOR_ELT(results, 4, wfilter);
  
  UNPROTECT(2);
  return results;
}

/* This function calculates the matrix of PHI[Jk](x[i]) for every k value where
 * phi[Jk] is non-null. The index i corresponds to the (i+1)-th line of the matrix.
 */
SEXP C_PHImat(SEXP x, SEXP J, SEXP family, SEXP taps, SEXP prec, SEXP periodic){
  
  int i, j, kdiff1, kmax, kmin, lkmin, lkmax, n, N, p, rJ, rper, rprec;
  double rphisl, rphisr, px, x1 = NA_REAL, xn = NA_REAL;
  double *rphimat1, *phi, *rwfilter, *rphimat2, *rx, *prod, *tmp;
  SEXP phimat;
  
  n = length(x);
  N = INTEGER(taps)[0];
  
  rJ = INTEGER(J)[0];
  rper = INTEGER(periodic)[0];
  rprec = INTEGER(prec)[0];
  rx = REAL(x);
  p = pow(2, rJ);
  
  // ----- Defining {min(x) --> x1} and {max(x) --> xn} ----- //
  Range(&x1, &xn, rx, n);
  
  if(ISNA(x1))
    error("Check your data. The observations should be real valued.");
  // ----- min(x) and max(x) defined ----- //
  
  SEXP wutils = PROTECT(WavUtilities(family, taps));
  SEXP phisl = VECTOR_ELT(wutils, 0);
  SEXP phisr = VECTOR_ELT(wutils, 1);
  SEXP wfilter = VECTOR_ELT(wutils, 4);
  rphisl = INTEGER(phisl)[0];
  rphisr = INTEGER(phisr)[0];
  rwfilter = REAL(wfilter);
  
  kmax = floor(p*xn - rphisl);
  kmin = ceil(p*x1 - rphisr + 1e-9);
  kdiff1 = kmax - kmin + 1;
  
  if(rper){
    rphimat1 = (double *) R_alloc(n*kdiff1, sizeof(double));
    PROTECT(phimat = allocMatrix(REALSXP, n, p));
    rphimat2 = REAL(phimat);
  }
  else{
    PROTECT(phimat = allocMatrix(REALSXP, n, kdiff1));
    rphimat1 = REAL(phimat);
  }
  
  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  
  // Let's fill the matrix PHI[Jk](x[i])!
  for(i = 0; i < n; i++){
    
    if(!R_FINITE(rx[i])){
      warning("At least one observation is not finite and its associated line in the PHI matrix will be set as NA.");
      for(j = 0; j < kdiff1; j++)
        rphimat1[i + n*j] = NA_REAL;
      continue;
    }
    
    px = p*rx[i];
    // ----- Calculating phi[0k](px) ----- //
    PhiVec(phi, px, rwfilter, N, rprec, prod, tmp);
    // ----- phi[0k](px) calculated --> phi ----- //
    
    lkmax = floor(px - rphisl);
    lkmin = lkmax - N + 2;
    
    // ----- Putting the phi[Jk](x[i]) in the i-th row of the matrix rphimat1 ----- //
    if(kmin == lkmin){
      for(j = 0; j < (N - 1); j++){
        rphimat1[i + n*j] = sqrt(p) * phi[j];
      }
      for(j = (N - 1); j < kdiff1; j++){
        rphimat1[i + n*j] = 0.0;
      }
    }
    else{
      for(j = 0; j < (lkmin - kmin); j++)
        rphimat1[i + n*j] = 0.0;
      
      for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
        rphimat1[i + n*j] = sqrt(p) * phi[j + kmin - lkmin];
      
      if(kmax > lkmax){
        for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
          rphimat1[i + n*j] = 0.0;
      }
    }
    // ----- phi[Jk](x[i]) put in the i-th row of the matrix rphimat1 ----- //
  }
  
  if(rper)
    Periodize(rphimat2, rphimat1, n, kmin, kmax, p);
  
  UNPROTECT(2);
  return phimat;
}

/* This function calculates the matrix of PSI[Jk](x[i]) for every k value where
 * psi[Jk] is non-null. The index i corresponds to the (i+1)-th line of the matrix.
 */
SEXP C_PSImat(SEXP x, SEXP J, SEXP family, SEXP taps, SEXP prec, SEXP periodic){
  
  int i, j, kdiff1, kmax, kmin, lkmin, lkmax, n, N, p, rJ, rper, rprec;
  double rpsisl, rpsisr, px, x1 = NA_REAL, xn = NA_REAL;
  double *psi, *rpsimat1, *rpsimat2, *rwfilter, *rx, *prod, *tmp;
  SEXP psimat;
  
  n = length(x);
  N = INTEGER(taps)[0];
  
  rJ = INTEGER(J)[0];
  rper = INTEGER(periodic)[0];
  rprec = INTEGER(prec)[0];
  rx = REAL(x);
  
  // ----- Defining {min(x) --> x1} and {max(x) --> xn} ----- //
  Range(&x1, &xn, rx, n);
  
  if(ISNA(x1))
    error("Check your data. The observations should be real valued.");
  // ----- min(x) and max(x) defined ----- //
  
  SEXP wutils = PROTECT(WavUtilities(family, taps));
  SEXP psisl = VECTOR_ELT(wutils, 2);
  SEXP psisr = VECTOR_ELT(wutils, 3);
  SEXP wfilter = VECTOR_ELT(wutils, 4);
  rpsisl = INTEGER(psisl)[0];
  rpsisr = INTEGER(psisr)[0];
  rwfilter = REAL(wfilter);
  
  p = pow(2, rJ);
  kmax = floor(p*xn - rpsisl);
  kmin = ceil(p*x1 - rpsisr + 1e-9);
  kdiff1 = kmax - kmin + 1;
  
  psi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  
  if(rper){
    rpsimat1 = (double *) R_alloc((n*kdiff1), sizeof(double));
    PROTECT(psimat = allocMatrix(REALSXP, n, p));
    rpsimat2 = REAL(psimat);
  }
  else{
    PROTECT(psimat = allocMatrix(REALSXP, n, kdiff1));
    rpsimat1 = REAL(psimat);
  }
  
  // Let's fill the matrix PSI[Jk](x[i])!
  for(i = 0; i < n; i++){
    
    if(!R_FINITE(rx[i])){
      warning("At least one observation is not finite and its associated line in the PSI matrix will be set as NA.");
      for(j = 0; j < kdiff1; j++)
        rpsimat1[i + n*j] = NA_REAL;
      continue;
    }
    
    px = p * rx[i];
    lkmax = floor(px - rpsisl);
    lkmin = lkmax - N + 2;
    
    // ----- Calculating phi[0k](px) ----- //
    PsiVec(psi, px, rwfilter, N, rprec, lkmin, prod, tmp);
    // ----- psi[0k](px) calculated --> psi ----- //
    
    // ----- Putting the psi[Jk](x[i]) in the i-th row of the matrix rpsimat1 ----- //
    if(kmin == lkmin){
      for(j = 0; j < (N - 1); j++)
        rpsimat1[i + n*j] = sqrt(p) * psi[j];
      
      for(j = (N - 1); j < kdiff1; j++)
        rpsimat1[i + n*j] = 0.0;
    }
    else{
      for(j = 0; j < (lkmin - kmin); j++)
        rpsimat1[i + n*j] = 0.0;
      
      for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
        rpsimat1[i + n*j] = sqrt(p) * psi[j + kmin - lkmin];
      
      if(kmax > lkmax){
        for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
          rpsimat1[i + n*j] = 0.0;
      }
    }
    // ----- psi[Jk](x[i]) put in the i-th row of the matrix rpsimat1 ----- //
  }
  
  if(rper)
    Periodize(rpsimat2, rpsimat1, n, kmin, kmax, p);
  
  UNPROTECT(2);
  return psimat;
}

SEXP PHImat(double *x, int n, int p, double *filter, int N, int prec, int kmin, int kmax, int phisl, int phisr, int periodic){
  
  int i, j, kdiff1, lkmin, lkmax;
  double px, *rphimat1, *phi, *rphimat2, *prod, *tmp;
  SEXP phimat;
  
  //n = length(x);
  //N = length(filter);
  
  //rfilter = REAL(filter);
  //rper = INTEGER(periodic)[0];
  //rphilh = INTEGER(philh)[0];
  //rphirh = INTEGER(phirh)[0];
  //rprec = INTEGER(prec)[0];
  //rx = REAL(x);
  
  kdiff1 = kmax - kmin + 1;
    
  if(periodic){
    rphimat1 = (double *) R_alloc(n*kdiff1, sizeof(double));
    PROTECT(phimat = allocMatrix(REALSXP, n, p));
    rphimat2 = REAL(phimat);
  }
  else{
    PROTECT(phimat = allocMatrix(REALSXP, n, kdiff1));
    rphimat1 = REAL(phimat);
  }
  
  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  
  // Let's fill the matrix PHI[Jk](x[i])!
  for(i = 0; i < n; i++){
    
    if(!R_FINITE(x[i])){
      for(j = 0; j < kdiff1; j++)
        rphimat1[i + n*j] = NA_REAL;
      continue;
    }
    
    px = p*x[i];
    // ----- Calculating phi[0k](px) ----- //
    PhiVec(phi, px, filter, N, prec, prod, tmp);
    // ----- phi[0k](px) calculated --> phi ----- //
    
    lkmax = floor(px - phisl);
    lkmin = lkmax - N + 2;
    
    // ----- Putting the phi[Jk](x[i]) in the i-th row of the matrix rphimat1 ----- //
    if(kmin == lkmin){
      for(j = 0; j < (N - 1); j++){
        rphimat1[i + n*j] = sqrt(p) * phi[j];
      }
      for(j = (N - 1); j < kdiff1; j++){
        rphimat1[i + n*j] = 0.0;
      }
    }
    else{
      for(j = 0; j < (lkmin - kmin); j++)
        rphimat1[i + n*j] = 0.0;
      
      for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
        rphimat1[i + n*j] = sqrt(p) * phi[j + kmin - lkmin];
      
      if(kmax > lkmax){
        for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
          rphimat1[i + n*j] = 0.0;
      }
    }
    // ----- phi[Jk](x[i]) put in the i-th row of the matrix rphimat1 ----- //
  }
  
  if(periodic)
    Periodize(rphimat2, rphimat1, n, kmin, kmax, p);
  
  UNPROTECT(1);
  return phimat;
}

SEXP PSImat(double *x, int n, int p, double *filter, int N, int prec, int kmin, int kmax, int psilh, int psirh, int periodic){
    
    int i, j, kdiff1, lkmin, lkmax;
    double px, *psi, *rpsimat1, *rpsimat2, *prod, *tmp;
    SEXP psimat;
    
    //n = length(x);
    //N = length(filter);
    
    //rfilter = REAL(filter);
    //rkmax = INTEGER(kmax)[0];
    //rkmin = INTEGER(kmin)[0];
    //rp = INTEGER(p)[0];
    //rper = INTEGER(periodic)[0];
    //rpsilh = INTEGER(psilh)[0];
    //rpsirh = INTEGER(psirh)[0];
    //rprec = INTEGER(prec)[0];
    //rx = REAL(x);
    
    kdiff1 = kmax - kmin + 1;
    
    psi = (double *) R_alloc((N - 1), sizeof(double));
    prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
    tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
    
    if(periodic){
        rpsimat1 = (double *) R_alloc(n*kdiff1, sizeof(double));
        PROTECT(psimat = allocMatrix(REALSXP, n, p));
        rpsimat2 = REAL(psimat);
    }
    else{
        PROTECT(psimat = allocMatrix(REALSXP, n, kdiff1));
        rpsimat1 = REAL(psimat);
    }
    
    for(i = 0; i < n; i++){
        
        if(!R_FINITE(x[i])){
            for(j = 0; j < kdiff1; j++)
                rpsimat1[i + n*j] = NA_REAL;
            continue;
        }
        
        px = p * x[i];
        lkmax = floor(px - psilh);
        lkmin = lkmax - N + 2;
        
        // ----- Calculating phi[0k](px) ----- //
        PsiVec(psi, px, filter, N, prec, lkmin, prod, tmp);
        // ----- psi[0k](px) calculated --> psi ----- //
        
        // ----- Putting the psi[Jk](x[i]) in the i-th row of the matrix rpsimat1 ----- //
        if(kmin == lkmin){
            for(j = 0; j < (N - 1); j++)
                rpsimat1[i + n*j] = sqrt(p) * psi[j];
            
            for(j = (N - 1); j < kdiff1; j++)
                rpsimat1[i + n*j] = 0.0;
        }
        else{
            for(j = 0; j < (lkmin - kmin); j++)
                rpsimat1[i + n*j] = 0.0;
            
            for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
                rpsimat1[i + n*j] = sqrt(p) * psi[j + kmin - lkmin];
            
            if(kmax > lkmax){
                for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
                    rpsimat1[i + n*j] = 0.0;
            }
        }
        // ----- psi[Jk](x[i]) put in the i-th row of the matrix rpsimat1 ----- //
    }
    
    if(periodic)
        Periodize(rpsimat2, rpsimat1, n, kmin, kmax, p);
    
    UNPROTECT(1);
    return psimat;
}

SEXP C_WavBasis(SEXP x, SEXP J0, SEXP J, SEXP family, SEXP taps, SEXP prec, SEXP periodic){
  
  double *rx, *rwfilt, x1 = NA_REAL, xn = NA_REAL;
  int i, n, N, rJ0, rJ, kmax, kmin, p, rper, rphisl, rphisr, rpsisl, rpsisr, rprec;
  
  rJ0 = INTEGER(J0)[0];
  rJ = INTEGER(J)[0];
  rprec = INTEGER(prec)[0];
  rper = INTEGER(periodic)[0];
  
  if(rJ0 > rJ)
    error("The coarsest level can't be greater than the finest level.");
  
  p = pow(2, rJ0);
  
  rx = REAL(x);
  n = LENGTH(x);
  N = INTEGER(taps)[0];
  
  // ----- Defining {min(x) --> x1} and {max(x) --> xn} ----- //
  Range(&x1, &xn, rx, n);
  
  if(ISNA(x1))
    error("Check your data. The observations should be real valued.");
  // ----- min(x) and max(x) defined ----- //
  
  SEXP wutils = PROTECT(WavUtilities(family, taps));
  SEXP phisl = VECTOR_ELT(wutils, 0);
  SEXP phisr = VECTOR_ELT(wutils, 1);
  SEXP psisl = VECTOR_ELT(wutils, 2);
  SEXP psisr = VECTOR_ELT(wutils, 3);
  SEXP wfilt = VECTOR_ELT(wutils, 4);
  rphisl = INTEGER(phisl)[0];
  rphisr = INTEGER(phisr)[0];
  rpsisl = INTEGER(psisl)[0];
  rpsisr = INTEGER(psisr)[0];
  rwfilt = REAL(wfilt);
  
  kmax = floor(p*xn - rphisl);
  kmin = ceil(p*x1 - rphisr + 1e-9);
  
  SEXP wmat, pmat;
  
  if(rJ0 == rJ){
    
    PROTECT(pmat = PHImat(rx, n, p, rwfilt, N, rprec, kmin, kmax, rphisl, rphisr, rper));
    wmat = pmat;
    
  }
  else{
    
    double *rpmat, *rwmat;
    int j, k, nc0, nc1, ncw;
    
    if(rper)
      ncw = pow(2, rJ);
    else{
      ncw = (kmax - kmin + 1);
      for(k = rJ0; k < rJ; k++){
        p = pow(2, k);
        kmax = floor(p*xn - rpsisl);
        kmin = ceil(p*x1 - rpsisr + 1e-9);
        ncw += (kmax - kmin + 1);
      }
    }
    
    PROTECT(wmat = allocMatrix(REALSXP, n, ncw));
    rwmat = REAL(wmat);
    
    p = pow(2, rJ0);
    kmax = floor(p*xn - rphisl);
    kmin = ceil(p*x1 - rphisr + 1e-9);
    
    PROTECT(pmat = PHImat(rx, n, p, rwfilt, N, rprec, kmin, kmax, rphisl, rphisr, rper));
    rpmat = REAL(pmat);
    
    nc0 = 0;
    nc1 = rper ? p : (kmax - kmin + 1);
    
    for(i = 0; i < n; i++)
      for(j = nc0; j < nc1; j++)
        rwmat[i + n*j] = rpmat[i + n*(j - nc0)];
    
    UNPROTECT(1);
    
    for(k = rJ0; k < rJ; k++){
      p = pow(2, k);
      kmax = floor(p*xn - rpsisl);
      kmin = ceil(p*x1 - rpsisr + 1e-9);
      PROTECT(pmat = PSImat(rx, n, p, rwfilt, N, rprec, kmin, kmax, rpsisl, rpsisr, rper));
      rpmat = REAL(pmat);
      
      nc0 = nc1;
      nc1 += (rper ? p : (kmax - kmin + 1));
      
      for(i = 0; i < n; i++)
        for(j = nc0; j < nc1; j++)
          rwmat[i + n*j] = rpmat[i + n*(j - nc0)];
          
      UNPROTECT(1);
      
    }
    
  }
  
  UNPROTECT(2);
  return wmat;
  
}
