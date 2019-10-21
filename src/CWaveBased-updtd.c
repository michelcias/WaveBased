
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
    w = w + 1e-09;
  else if(w == 1)
    w = w - 1rwf9;
  
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
 * pywavelets website (http://wavelets.pybytes.com).
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
      rwfilter[0]  =  0.11154074335008017;
      rwfilter[1]  =  0.4946238903983854;
      rwfilter[2]  =  0.7511339080215775;
      rwfilter[3]  =  0.3152503517092432;
      rwfilter[4]  = -0.22626469396516913;
      rwfilter[5]  = -0.12976686756709563;
      rwfilter[6]  =  0.09750160558707936;
      rwfilter[7]  =  0.02752286553001629;
      rwfilter[8]  = -0.031582039318031156;
      rwfilter[9]  =  0.0005538422009938016;
      rwfilter[10] =  0.004777257511010651;
      rwfilter[11] = -0.00107730108499558;
    }
    else if(rtaps == 14){
      rwfilter[0]  =  0.07785205408506236;
      rwfilter[1]  =  0.39653931948230575;
      rwfilter[2]  =  0.7291320908465551;
      rwfilter[3]  =  0.4697822874053586;
      rwfilter[4]  = -0.14390600392910627;
      rwfilter[5]  = -0.22403618499416572;
      rwfilter[6]  =  0.07130921926705004;
      rwfilter[7]  =  0.0806126091510659;
      rwfilter[8]  = -0.03802993693503463;
      rwfilter[9]  = -0.01657454163101562;
      rwfilter[10] =  0.012550998556013784;
      rwfilter[11] =  0.00042957797300470274;
      rwfilter[12] = -0.0018016407039998328;
      rwfilter[13] =  0.0003537138000010399;
    }
    else if(rtaps == 16){
      rwfilter[0]  =  0.05441584224308161;
      rwfilter[1]  =  0.3128715909144659;
      rwfilter[2]  =  0.6756307362980128;
      rwfilter[3]  =  0.5853546836548691;
      rwfilter[4]  = -0.015829105256023893;
      rwfilter[5]  = -0.2840155429624281;
      rwfilter[6]  =  0.00047248457399797254;
      rwfilter[7]  =  0.128747426620186;
      rwfilter[8]  = -0.01736930100202211;
      rwfilter[9]  = -0.04408825393106472;
      rwfilter[10] =  0.013981027917015516;
      rwfilter[11] =  0.008746094047015655;
      rwfilter[12] = -0.00487035299301066;
      rwfilter[13] = -0.0003917403729959771;
      rwfilter[14] =  0.0006754494059985568;
      rwfilter[15] = -0.00011747678400228192;
    }
    else if(rtaps == 18){
      rwfilter[0]  =  0.03807794736316728;
      rwfilter[1]  =  0.24383467463766728;
      rwfilter[2]  =  0.6048231236767786;
      rwfilter[3]  =  0.6572880780366389;
      rwfilter[4]  =  0.13319738582208895;
      rwfilter[5]  = -0.29327378327258685;
      rwfilter[6]  = -0.09684078322087904;
      rwfilter[7]  =  0.14854074933476008;
      rwfilter[8]  =  0.030725681478322865;
      rwfilter[9]  = -0.06763282905952399;
      rwfilter[10] =  0.00025094711499193845;
      rwfilter[11] =  0.022361662123515244;
      rwfilter[12] = -0.004723204757894831;
      rwfilter[13] = -0.004281503681904723;
      rwfilter[14] =  0.0018476468829611268;
      rwfilter[15] =  0.00023038576399541288;
      rwfilter[16] = -0.0002519631889981789;
      rwfilter[17] =  3.9347319995026124e-05;
    }
    else if(rtaps == 20){
      rwfilter[0]  =  0.026670057900950818;
      rwfilter[1]  =  0.18817680007762133;
      rwfilter[2]  =  0.5272011889309198;
      rwfilter[3]  =  0.6884590394525921;
      rwfilter[4]  =  0.2811723436604265;
      rwfilter[5]  = -0.24984642432648865;
      rwfilter[6]  = -0.19594627437659665;
      rwfilter[7]  =  0.12736934033574265;
      rwfilter[8]  =  0.09305736460380659;
      rwfilter[9]  = -0.07139414716586077;
      rwfilter[10] = -0.02945753682194567;
      rwfilter[11] =  0.03321267405893324;
      rwfilter[12] =  0.0036065535669883944;
      rwfilter[13] = -0.010733175482979604;
      rwfilter[14] =  0.0013953517469940798;
      rwfilter[15] =  0.00199240529499085;
      rwfilter[16] = -0.0006858566950046825;
      rwfilter[17] = -0.0001164668549943862;
      rwfilter[18] =  9.358867000108985e-05;
      rwfilter[19] = -1.326420300235487e-05;
    }
    else if(rtaps == 22){
      rwfilter[0]  =  0.01869429776147044;
      rwfilter[1]  =  0.1440670211506196;
      rwfilter[2]  =  0.44989976435603013;
      rwfilter[3]  =  0.6856867749161785;
      rwfilter[4]  =  0.41196436894789695;
      rwfilter[5]  = -0.16227524502747828;
      rwfilter[6]  = -0.27423084681792875;
      rwfilter[7]  =  0.06604358819669089;
      rwfilter[8]  =  0.14981201246638268;
      rwfilter[9]  = -0.04647995511667613;
      rwfilter[10] = -0.06643878569502022;
      rwfilter[11] =  0.03133509021904531;
      rwfilter[12] =  0.02084090436018004;
      rwfilter[13] = -0.015364820906201324;
      rwfilter[14] = -0.0033408588730145018;
      rwfilter[15] =  0.004928417656058778;
      rwfilter[16] = -0.00030859285881515924;
      rwfilter[17] = -0.0008930232506662366;
      rwfilter[18] =  0.00024915252355281426;
      rwfilter[19] =  5.443907469936638e-05;
      rwfilter[10] = -3.463498418698379e-05;
      rwfilter[21] =  4.494274277236352e-06;
    }
    else if(rtaps == 24){
      rwfilter[0] =  0.013112257957229239;
      rwfilter[1] =  0.10956627282118277;
      rwfilter[2] =  0.3773551352142041;
      rwfilter[3] =  0.6571987225792911;
      rwfilter[4] =  0.5158864784278007;
      rwfilter[5] = -0.04476388565377762;
      rwfilter[6] = -0.31617845375277914;
      rwfilter[7] = -0.023779257256064865;
      rwfilter[8] =  0.18247860592758275;
      rwfilter[9] =  0.0053595696743599965;
      rwfilter[10] = -0.09643212009649671;
      rwfilter[11] =  0.010849130255828966;
      rwfilter[12] =  0.04154627749508764;
      rwfilter[13] = -0.01221864906974642;
      rwfilter[14] = -0.012840825198299882;
      rwfilter[15] =  0.006711499008795549;
      rwfilter[16] =  0.0022486072409952287;
      rwfilter[17] = -0.0021795036186277044;
      rwfilter[18] =  6.5451282125215034e-06;
      rwfilter[19] =  0.0003886530628209267;
      rwfilter[20] = -8.850410920820318e-05;
      rwfilter[21] = -2.4241545757030318e-05;
      rwfilter[22] =  1.2776952219379579e-05;
      rwfilter[23] = -1.5290717580684923e-06;
    }
    else if(rtaps == 26){
      rwfilter[0]  =  0.009202133538962279;
      rwfilter[1]  =  0.08286124387290195;
      rwfilter[2]  =  0.3119963221604349;
      rwfilter[3]  =  0.6110558511587811;
      rwfilter[4]  =  0.5888895704312119;
      rwfilter[5]  =  0.086985726179645;
      rwfilter[6]  = -0.31497290771138414;
      rwfilter[7]  = -0.12457673075080665;
      rwfilter[8]  =  0.17947607942935084;
      rwfilter[9]  =  0.07294893365678874;
      rwfilter[10] = -0.10580761818792761;
      rwfilter[11] = -0.026488406475345658;
      rwfilter[12] =  0.056139477100276156;
      rwfilter[13] =  0.002379972254052227;
      rwfilter[14] = -0.02383142071032781;
      rwfilter[15] =  0.003923941448795577;
      rwfilter[16] =  0.007255589401617119;
      rwfilter[17] = -0.002761911234656831;
      rwfilter[18] = -0.0013156739118922766;
      rwfilter[19] =  0.000932326130867249;
      rwfilter[20] =  4.9251525126285676e-05;
      rwfilter[21] = -0.0001651289885565057;
      rwfilter[22] =  3.067853757932436e-05;
      rwfilter[23] =  1.0441930571407941e-05;
      rwfilter[24] = -4.700416479360808e-06;
      rwfilter[25] =  5.2200350984548e-07;
    }
    else if(rtaps == 28){
      rwfilter[0]  =  0.0064611534600864905;
      rwfilter[1]  =  0.062364758849384874;
      rwfilter[2]  =  0.25485026779256437;
      rwfilter[3]  =  0.5543056179407709;
      rwfilter[4]  =  0.6311878491047198;
      rwfilter[5]  =  0.21867068775886594;
      rwfilter[6]  = -0.27168855227867705;
      rwfilter[7]  = -0.2180335299932165;
      rwfilter[8]  =  0.13839521386479153;
      rwfilter[9]  =  0.13998901658445695;
      rwfilter[10] = -0.0867484115681106;
      rwfilter[11] = -0.0715489555039835;
      rwfilter[12] =  0.05523712625925082;
      rwfilter[13] =  0.02698140830794797;
      rwfilter[14] = -0.030185351540353976;
      rwfilter[15] = -0.0056150495303375755;
      rwfilter[16] =  0.01278949326634007;
      rwfilter[17] = -0.0007462189892638753;
      rwfilter[18] = -0.003849638868019787;
      rwfilter[19] =  0.001061691085606874;
      rwfilter[20] =  0.0007080211542354048;
      rwfilter[21] = -0.00038683194731287514;
      rwfilter[22] = -4.177724577037067e-05;
      rwfilter[23] =  6.875504252695734e-05;
      rwfilter[24] = -1.0337209184568496e-05;
      rwfilter[25] = -4.389704901780418e-06;
      rwfilter[26] =  1.7249946753674012e-06;
      rwfilter[27] = -1.7871399683109222e-07;
    }
    else if(rtaps == 30){
      rwfilter[0]  =  0.004538537361577376;
      rwfilter[1]  =  0.04674339489275062;
      rwfilter[2]  =  0.20602386398692688;
      rwfilter[3]  =  0.4926317717079753;
      rwfilter[4]  =  0.6458131403572103;
      rwfilter[5]  =  0.33900253545462167;
      rwfilter[6]  = -0.19320413960907623;
      rwfilter[7]  = -0.28888259656686216;
      rwfilter[8]  =  0.06528295284876569;
      rwfilter[9]  =  0.19014671400708816;
      rwfilter[10] = -0.0396661765557336;
      rwfilter[11] = -0.11112093603713753;
      rwfilter[12] =  0.033877143923563204;
      rwfilter[13] =  0.054780550584559995;
      rwfilter[14] = -0.02576700732836694;
      rwfilter[15] = -0.020810050169636805;
      rwfilter[16] =  0.015083918027862582;
      rwfilter[17] =  0.005101000360422873;
      rwfilter[18] = -0.0064877345603061454;
      rwfilter[19] = -0.00024175649075894543;
      rwfilter[20] =  0.0019433239803823459;
      rwfilter[21] = -0.0003734823541372647;
      rwfilter[22] = -0.00035956524436229364;
      rwfilter[23] =  0.00015589648992055726;
      rwfilter[24] =  2.579269915531323e-05;
      rwfilter[25] = -2.8133296266037558e-05;
      rwfilter[26] =  3.3629871817363823e-06;
      rwfilter[27] =  1.8112704079399406e-06;
      rwfilter[28] = -6.316882325879451e-07;
      rwfilter[29] =  6.133359913303714e-08;
    }
    else if(rtaps == 32){
      rwfilter[0]  =  0.0031892209253436892;
      rwfilter[1]  =  0.03490771432362905;
      rwfilter[2]  =  0.1650642834886438;
      rwfilter[3]  =  0.43031272284545874;
      rwfilter[4]  =  0.6373563320829833;
      rwfilter[5]  =  0.44029025688580486;
      rwfilter[6]  = -0.08975108940236352;
      rwfilter[7]  = -0.3270633105274758;
      rwfilter[8]  = -0.02791820813292813;
      rwfilter[9]  =  0.21119069394696974;
      rwfilter[10] =  0.027340263752899923;
      rwfilter[11] = -0.13238830556335474;
      rwfilter[12] = -0.006239722752156254;
      rwfilter[13] =  0.07592423604445779;
      rwfilter[14] = -0.007588974368642594;
      rwfilter[15] = -0.036888397691556774;
      rwfilter[16] =  0.010297659641009963;
      rwfilter[17] =  0.013993768859843242;
      rwfilter[18] = -0.006990014563390751;
      rwfilter[19] = -0.0036442796214883506;
      rwfilter[20] =  0.00312802338120381;
      rwfilter[21] =  0.00040789698084934395;
      rwfilter[22] = -0.0009410217493585433;
      rwfilter[23] =  0.00011424152003843815;
      rwfilter[24] =  0.00017478724522506327;
      rwfilter[25] = -6.103596621404321e-05;
      rwfilter[26] = -1.394566898819319e-05;
      rwfilter[27] =  1.133660866126152e-05;
      rwfilter[28] = -1.0435713423102517e-06;
      rwfilter[29] = -7.363656785441815e-07;
      rwfilter[30] =  2.3087840868545578e-07;
      rwfilter[31] = -2.1093396300980412e-08;
    }
    else if(rtaps == 34){
      rwfilter[0]  =  0.00224180700103879;
      rwfilter[1]  =  0.025985393703623173;
      rwfilter[2]  =  0.13121490330791097;
      rwfilter[3]  =  0.3703507241528858;
      rwfilter[4]  =  0.6109966156850273;
      rwfilter[5]  =  0.5183157640572823;
      rwfilter[6]  =  0.027314970403312946;
      rwfilter[7]  = -0.32832074836418546;
      rwfilter[8]  = -0.12659975221599248;
      rwfilter[9]  =  0.19731058956508457;
      rwfilter[10] =  0.10113548917744287;
      rwfilter[11] = -0.12681569177849797;
      rwfilter[12] = -0.05709141963185808;
      rwfilter[13] =  0.08110598665408082;
      rwfilter[14] =  0.022312336178011833;
      rwfilter[15] = -0.04692243838937891;
      rwfilter[16] = -0.0032709555358783646;
      rwfilter[17] =  0.022733676583919053;
      rwfilter[18] = -0.0030429899813869555;
      rwfilter[19] = -0.008602921520347815;
      rwfilter[20] =  0.002967996691518064;
      rwfilter[21] =  0.0023012052421511474;
      rwfilter[22] = -0.001436845304805;
      rwfilter[23] = -0.00032813251941022427;
      rwfilter[24] =  0.0004394654277689454;
      rwfilter[25] = -2.5610109566546042e-05;
      rwfilter[26] = -8.204803202458212e-05;
      rwfilter[27] =  2.318681379876164e-05;
      rwfilter[28] =  6.990600985081294e-06;
      rwfilter[29] = -4.505942477225963e-06;
      rwfilter[30] =  3.0165496099963414e-07;
      rwfilter[31] =  2.9577009333187617e-07;
      rwfilter[32] = -8.423948446008154e-08;
      rwfilter[33] =  7.26749296856637e-09;
    }
    else if(rtaps == 36){
      rwfilter[0]  =  0.0015763102184365595;
      rwfilter[1]  =  0.01928853172409497;
      rwfilter[2]  =  0.10358846582214751;
      rwfilter[3]  =  0.31467894133619284;
      rwfilter[4]  =  0.5718268077650818;
      rwfilter[5]  =  0.571801654887122;
      rwfilter[6]  =  0.14722311196952223;
      rwfilter[7]  = -0.2936540407357981;
      rwfilter[8]  = -0.21648093400458224;
      rwfilter[9]  =  0.14953397556500755;
      rwfilter[10] =  0.16708131276294505;
      rwfilter[11] = -0.09233188415030412;
      rwfilter[12] = -0.10675224665906288;
      rwfilter[13] =  0.0648872162123582;
      rwfilter[14] =  0.05705124773905827;
      rwfilter[15] = -0.04452614190225633;
      rwfilter[16] = -0.023733210395336858;
      rwfilter[17] =  0.026670705926689853;
      rwfilter[18] =  0.006262167954438661;
      rwfilter[19] = -0.013051480946517112;
      rwfilter[20] =  0.00011863003387493042;
      rwfilter[21] =  0.004943343605456594;
      rwfilter[22] = -0.0011187326669886426;
      rwfilter[23] = -0.0013405962983313922;
      rwfilter[24] =  0.0006284656829644715;
      rwfilter[25] =  0.0002135815619103188;
      rwfilter[26] = -0.00019864855231101547;
      rwfilter[27] = -1.535917123021341e-07;
      rwfilter[28] =  3.741237880730847e-05;
      rwfilter[29] = -8.520602537423464e-06;
      rwfilter[30] = -3.3326344788769603e-06;
      rwfilter[31] =  1.768712983622886e-06;
      rwfilter[32] = -7.691632689865049e-08;
      rwfilter[33] = -1.1760987670250871e-07;
      rwfilter[34] =  3.06883586303703e-08;
      rwfilter[35] = -2.507934454941929e-09;
    }
    else if(rtaps == 38){
      rwfilter[0]  =  0.0011086697631864314;
      rwfilter[1]  =  0.01428109845082521;
      rwfilter[2]  =  0.08127811326580564;
      rwfilter[3]  =  0.26438843174202237;
      rwfilter[4]  =  0.5244363774668862;
      rwfilter[5]  =  0.6017045491300916;
      rwfilter[6]  =  0.2608949526521201;
      rwfilter[7]  = -0.22809139421653665;
      rwfilter[8]  = -0.28583863175723145;
      rwfilter[9]  =  0.07465226970806647;
      rwfilter[10] =  0.21234974330662043;
      rwfilter[11] = -0.03351854190320226;
      rwfilter[12] = -0.14278569504021468;
      rwfilter[13] =  0.02758435062488713;
      rwfilter[14] =  0.0869067555554507;
      rwfilter[15] = -0.026501236250778635;
      rwfilter[16] = -0.04567422627778492;
      rwfilter[17] =  0.021623767409452484;
      rwfilter[18] =  0.019375549889114482;
      rwfilter[19] = -0.013988388678695632;
      rwfilter[20] = -0.005866922281112195;
      rwfilter[21] =  0.007040747367080495;
      rwfilter[22] =  0.0007689543592242488;
      rwfilter[23] = -0.002687551800734441;
      rwfilter[24] =  0.00034180865344939543;
      rwfilter[25] =  0.0007358025205041731;
      rwfilter[26] = -0.0002606761356811995;
      rwfilter[27] = -0.00012460079173506306;
      rwfilter[28] =  8.711270467250443e-05;
      rwfilter[29] =  5.105950487090694e-06;
      rwfilter[30] = -1.664017629722462e-05;
      rwfilter[31] =  3.0109643163099385e-06;
      rwfilter[32] =  1.531931476697877e-06;
      rwfilter[33] = -6.86275565779811e-07;
      rwfilter[34] =  1.447088298804088e-08;
      rwfilter[35] =  4.636937775802368e-08;
      rwfilter[36] = -1.1164020670405678e-08;
      rwfilter[37] =  8.666848839034483e-10;
    }
    else if(rtaps == 40){
      rwfilter[0]  =  0.0007799536136659112;
      rwfilter[1]  =  0.010549394624937735;
      rwfilter[2]  =  0.06342378045900529;
      rwfilter[3]  =  0.21994211355113222;
      rwfilter[4]  =  0.4726961853103315;
      rwfilter[5]  =  0.6104932389378558;
      rwfilter[6]  =  0.36150229873889705;
      rwfilter[7]  = -0.13921208801128787;
      rwfilter[8]  = -0.3267868004335376;
      rwfilter[9]  = -0.016727088308801888;
      rwfilter[10] =  0.22829105082013823;
      rwfilter[11] =  0.039850246458519104;
      rwfilter[12] = -0.1554587507060453;
      rwfilter[13] = -0.024716827337521424;
      rwfilter[14] =  0.10229171917513397;
      rwfilter[15] =  0.005632246857685454;
      rwfilter[16] = -0.061722899624668884;
      rwfilter[17] =  0.0058746818113949465;
      rwfilter[18] =  0.03229429953011916;
      rwfilter[19] = -0.008789324924555765;
      rwfilter[20] = -0.013810526137727442;
      rwfilter[21] =  0.0067216273018096935;
      rwfilter[22] =  0.00442054238676635;
      rwfilter[23] = -0.003581494259744107;
      rwfilter[24] = -0.0008315621728772474;
      rwfilter[25] =  0.0013925596193045254;
      rwfilter[26] = -5.349759844340453e-05;
      rwfilter[27] = -0.0003851047486990061;
      rwfilter[28] =  0.00010153288973669777;
      rwfilter[29] =  6.774280828373048e-05;
      rwfilter[30] = -3.710586183390615e-05;
      rwfilter[31] = -4.376143862182197e-06;
      rwfilter[32] =  7.241248287663791e-06;
      rwfilter[33] = -1.0119940100181473e-06;
      rwfilter[34] = -6.847079596993149e-07;
      rwfilter[35] =  2.633924226266962e-07;
      rwfilter[36] =  2.0143220235374613e-10;
      rwfilter[37] = -1.814843248297622e-08;
      rwfilter[38] =  4.05612705554717e-09;
      rwfilter[39] = -2.998836489615753e-10;
    }
    else
      error("'taps = %d' is not allowed for 'Daublets'. For this family, only 2, 4, 6, ..., 38 and 40 taps are available.", rtaps);
  }
  else if(rfam == 2){
      if(rtaps == 4){
        rwfilter[0] = -0.12940952255092145;
        rwfilter[1] =  0.22414386804185735;
        rwfilter[2] =  0.836516303737469;
        rwfilter[3] =  0.48296291314469025;
      }
      if(rtaps == 6){
        rwfilter[0] =  0.035226291882100656;
        rwfilter[1] = -0.08544127388224149;
        rwfilter[2] = -0.13501102001039084;
        rwfilter[3] =  0.4598775021193313;
        rwfilter[4] =  0.8068915093133388;
        rwfilter[5] =  0.3326705529509569;
      }
      if(rtaps == 8){
        rwfilter[0] = -0.07576571478927333;
        rwfilter[1] = -0.02963552764599851;
        rwfilter[2] =  0.49761866763201545;
        rwfilter[3] =  0.8037387518059161;
        rwfilter[4] =  0.29785779560527736;
        rwfilter[5] = -0.09921954357684722;
        rwfilter[6] = -0.012603967262037833;
        rwfilter[7] =  0.0322231006040427;
      }
    else if(rtaps == 10){
      rwfilter[0]  =  0.027333068345077982;
      rwfilter[1]  =  0.029519490925774643;
      rwfilter[2]  = -0.039134249302383094;
      rwfilter[3]  =  0.1993975339773936;
      rwfilter[4]  =  0.7234076904024206;
      rwfilter[5]  =  0.6339789634582119;
      rwfilter[6]  =  0.01660210576452232;
      rwfilter[7]  = -0.17532808990845047;
      rwfilter[8]  = -0.021101834024758855;
      rwfilter[9]  =  0.019538882735286728;
    }
    else if(rtaps == 12){
      rwfilter[0]  =  0.015404109327027373;
      rwfilter[1]  =  0.0034907120842174702;
      rwfilter[2]  = -0.11799011114819057;
      rwfilter[3]  = -0.048311742585633;
      rwfilter[4]  =  0.4910559419267466;
      rwfilter[5]  =  0.787641141030194;
      rwfilter[6]  =  0.3379294217276218;
      rwfilter[7]  = -0.07263752278646252;
      rwfilter[8]  = -0.021060292512300564;
      rwfilter[9]  =  0.04472490177066578;
      rwfilter[10] =  0.0017677118642428036;
      rwfilter[11] = -0.007800708325034148;
    }
    else if(rtaps == 14){
      rwfilter[0]  =  0.002681814568257878;
      rwfilter[1]  = -0.0010473848886829163;
      rwfilter[2]  = -0.01263630340325193;
      rwfilter[3]  =  0.03051551316596357;
      rwfilter[4]  =  0.0678926935013727;
      rwfilter[5]  = -0.049552834937127255;
      rwfilter[6]  =  0.017441255086855827;
      rwfilter[7]  =  0.5361019170917628;
      rwfilter[8]  =  0.767764317003164;
      rwfilter[9]  =  0.2886296317515146;
      rwfilter[10] = -0.14004724044296152;
      rwfilter[11] = -0.10780823770381774;
      rwfilter[12] =  0.004010244871533663;
      rwfilter[13] =  0.010268176708511255;
    }
    else if(rtaps == 16){
      rwfilter[0]  =  0.0018899503327594609;
      rwfilter[1]  = -0.0003029205147213668;
      rwfilter[2]  = -0.01495225833704823;
      rwfilter[3]  =  0.003808752013890615;
      rwfilter[4]  =  0.049137179673607506;
      rwfilter[5]  = -0.027219029917056003;
      rwfilter[6]  = -0.05194583810770904;
      rwfilter[7]  =  0.3644418948353314;
      rwfilter[8]  =  0.7771857517005235;
      rwfilter[9]  =  0.4813596512583722;
      rwfilter[10] = -0.061273359067658524;
      rwfilter[11] = -0.1432942383508097;
      rwfilter[12] =  0.007607487324917605;
      rwfilter[13] =  0.03169508781149298;
      rwfilter[14] = -0.0005421323317911481;
      rwfilter[15] = -0.0033824159510061256;
    }
    else if(rtaps == 18){
      rwfilter[0]  =  0.0010694900329086053;
      rwfilter[1]  = -0.0004731544986800831;
      rwfilter[2]  = -0.010264064027633142;
      rwfilter[3]  =  0.008859267493400484;
      rwfilter[4]  =  0.06207778930288603;
      rwfilter[5]  = -0.018233770779395985;
      rwfilter[6]  = -0.19155083129728512;
      rwfilter[7]  =  0.035272488035271894;
      rwfilter[8]  =  0.6173384491409358;
      rwfilter[9]  =  0.717897082764412;
      rwfilter[10] =  0.238760914607303;
      rwfilter[11] = -0.05456895843083407;
      rwfilter[12] =  0.0005834627461258068;
      rwfilter[13] =  0.03022487885827568;
      rwfilter[14] = -0.01152821020767923;
      rwfilter[15] = -0.013271967781817119;
      rwfilter[16] =  0.0006197808889855868;
      rwfilter[17] =  0.0014009155259146807;
    }
    else if(rtaps == 20){
      rwfilter[0]  =  0.0007701598091144901;
      rwfilter[1]  =  9.563267072289475e-05;
      rwfilter[2]  = -0.008641299277022422;
      rwfilter[3]  = -0.0014653825813050513;
      rwfilter[4]  =  0.0459272392310922;
      rwfilter[5]  =  0.011609893903711381;
      rwfilter[6]  = -0.15949427888491757;
      rwfilter[7]  = -0.07088053578324385;
      rwfilter[8]  =  0.47169066693843925;
      rwfilter[9]  =  0.7695100370211071;
      rwfilter[10] =  0.38382676106708546;
      rwfilter[11] = -0.03553674047381755;
      rwfilter[12] = -0.0319900568824278;
      rwfilter[13] =  0.04999497207737669;
      rwfilter[14] =  0.005764912033581909;
      rwfilter[15] = -0.02035493981231129;
      rwfilter[16] = -0.0008043589320165449;
      rwfilter[17] =  0.004593173585311828;
      rwfilter[18] =  5.7036083618494284e-05;
      rwfilter[19] = -0.0004593294210046588;
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
      rwfilter[0] = -0.0727326195128539;
      rwfilter[1] =  0.3378976624578092;
      rwfilter[2] =  0.8525720202122554;
      rwfilter[3] =  0.38486484686420286;
      rwfilter[4] = -0.0727326195128539;
      rwfilter[5] = -0.01565572813546454;
    }
    else if(rtaps == 12){
      rwfilter[0]  =  0.016387336463522112;
      rwfilter[1]  = -0.04146493678175915;
      rwfilter[2]  = -0.06737255472196302;
      rwfilter[3]  =  0.3861100668211622;
      rwfilter[4]  =  0.8127236354455423;
      rwfilter[5]  =  0.41700518442169254;
      rwfilter[6]  = -0.0764885990783064;
      rwfilter[7]  = -0.0594344186464569;
      rwfilter[8]  =  0.023680171946334084;
      rwfilter[9]  =  0.0056114348193944995;
      rwfilter[10] = -0.0018232088707029932;
      rwfilter[11] = -0.0007205494453645122;
    }
    else if(rtaps == 18){
      rwfilter[0]  = -0.003793512864491014;
      rwfilter[1]  =  0.007782596427325418;
      rwfilter[2]  =  0.023452696141836267;
      rwfilter[3]  = -0.0657719112818555;
      rwfilter[4]  = -0.06112339000267287;
      rwfilter[5]  =  0.4051769024096169;
      rwfilter[6]  =  0.7937772226256206;
      rwfilter[7]  =  0.42848347637761874;
      rwfilter[8]  = -0.07179982161931202;
      rwfilter[9]  = -0.08230192710688598;
      rwfilter[10] =  0.03455502757306163;
      rwfilter[11] =  0.015880544863615904;
      rwfilter[12] = -0.00900797613666158;
      rwfilter[13] = -0.0025745176887502236;
      rwfilter[14] =  0.0011175187708906016;
      rwfilter[15] =  0.0004662169601128863;
      rwfilter[16] = -7.098330313814125e-05;
      rwfilter[17] = -3.459977283621256e-05;
    }
    else if(rtaps == 24){
      rwfilter[0]  =  0.0008923136685823146;
      rwfilter[1]  = -0.0016294920126017326;
      rwfilter[2]  = -0.0073461663276420935;
      rwfilter[3]  =  0.016068943964776348;
      rwfilter[4]  =  0.026682300156053072;
      rwfilter[5]  = -0.08126669968087875;
      rwfilter[6]  = -0.05607731331675481;
      rwfilter[7]  =  0.41530840703043026;
      rwfilter[8]  =  0.782238930920499;
      rwfilter[9]  =  0.4343860564914685;
      rwfilter[10] = -0.06662747426342504;
      rwfilter[11] = -0.09622044203398798;
      rwfilter[12] =  0.03933442712333749;
      rwfilter[13] =  0.025082261844864097;
      rwfilter[14] = -0.015211731527946259;
      rwfilter[15] = -0.00565828668661072;
      rwfilter[16] =  0.003751436157278457;
      rwfilter[17] =  0.0012665619292989445;
      rwfilter[18] = -0.0005890207562443383;
      rwfilter[19] = -0.00025997455248771324;
      rwfilter[20] =  6.233903446100713e-05;
      rwfilter[21] =  3.1229875865345646e-05;
      rwfilter[22] = -3.2596802368833675e-06;
      rwfilter[23] = -1.7849850030882614e-06;
    }
    else if(rtaps == 30){
      rwfilter[0]  = -0.00021208083980379827;
      rwfilter[1]  =  0.00035858968789573785;
      rwfilter[2]  =  0.0021782363581090178;
      rwfilter[3]  = -0.004159358781386048;
      rwfilter[4]  = -0.010131117519849788;
      rwfilter[5]  =  0.023408156785839195;
      rwfilter[6]  =  0.02816802897093635;
      rwfilter[7]  = -0.09192001055969624;
      rwfilter[8]  = -0.05204316317624377;
      rwfilter[9]  =  0.4215662066908515;
      rwfilter[10] =  0.7742896036529562;
      rwfilter[11] =  0.4379916261718371;
      rwfilter[12] = -0.06203596396290357;
      rwfilter[13] = -0.10557420870333893;
      rwfilter[14] =  0.0412892087501817;
      rwfilter[15] =  0.03268357426711183;
      rwfilter[16] = -0.01976177894257264;
      rwfilter[17] = -0.009164231162481846;
      rwfilter[18] =  0.006764185448053083;
      rwfilter[19] =  0.0024333732126576722;
      rwfilter[20] = -0.0016628637020130838;
      rwfilter[21] = -0.0006381313430451114;
      rwfilter[22] =  0.00030225958181306315;
      rwfilter[23] =  0.00014054114970203437;
      rwfilter[24] = -4.134043227251251e-05;
      rwfilter[25] = -2.1315026809955787e-05;
      rwfilter[26] =  3.7346551751414047e-06;
      rwfilter[27] =  2.0637618513646814e-06;
      rwfilter[28] = -1.6744288576823017e-07;
      rwfilter[29] = -9.517657273819165e-08;
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
  kmin = ceil(p*x1 - rphisr + 1e-09);
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
  kmin = ceil(p*x1 - rpsisr + 1e-09);
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
  kmin = ceil(p*x1 - rphisr + 1e-09);
  
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
        kmin = ceil(p*x1 - rpsisr + 1e-09);
        ncw += (kmax - kmin + 1);
      }
    }
    
    PROTECT(wmat = allocMatrix(REALSXP, n, ncw));
    rwmat = REAL(wmat);
    
    p = pow(2, rJ0);
    kmax = floor(p*xn - rphisl);
    kmin = ceil(p*x1 - rphisr + 1e-09);
    
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
      kmin = ceil(p*x1 - rpsisr + 1e-09);
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
