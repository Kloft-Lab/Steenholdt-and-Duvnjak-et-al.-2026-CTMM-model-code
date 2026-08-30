;; 1. Based on: run60
;; 2. Description: final model, initials identical to final estimates (for sharing)
;; x1. Author: Zrinka Duvnjak
$PROB CTMM Model with Time-Varying Dosing Delay Hazard on Flaring

$INPUT ID GA=TIME PGA EVID GALASTD DOSINT SOURCE PREG PDV DV PPGA AMT CMT FIRSTOBS INITGA
$DATA ../../data/data_CTMM_new_data.csv IGNORE=@ IGNORE=(TIME.GT.70)

$SUBROUTINES ADVAN13 TOL=9

$MODEL
  COMP=(REMISSION) ; CMT=1
  COMP=(FLARE)     ; CMT=2


$PK
  ; Baseline rates (with log-normal random effect on flaring)
  L12_BASE = THETA(1) * EXP(ETA(1)) 
  L21_BASE = THETA(2) 

  ; Parameter for the effect of elapsed dosing intervals
  EFF_DEV = THETA(3)


$DES
  ; Calculate continuous Time Since Last Dose 
  ; 'T' is the continuous integration clock (which represents GA here)
  TSLD = T - GALASTD

  ; Safety catch: prevent negative time if solver overshoots backwards
  IF (TSLD .LT. 0) TSLD = 0

  ; Calculate continuous number of dosing intervals passed
  INT_PASSED = TSLD / DOSINT

  ; Calculate how many intervals late the patient is
  DELAYINT = INT_PASSED - 1

  ; Safety catch: If they are not late yet, delay is zero
  IF (DELAYINT .LT. 0) DELAYINT = 0

  ; Calculate the time-varying flaring rate
  L12 = L12_BASE * EXP(EFF_DEV * DELAYINT)
  L21 = L21_BASE 

  ; CTMC Differential Equations
  DADT(1) = -L12 * A(1) + L21 * A(2)
  DADT(2) =  L12 * A(1) - L21 * A(2)


$ERROR
  ; Pass the ODE compartment amounts into printable custom variables
  P_REM = A(1)
  P_FLARE = A(2)

  ; Create a temporary variable to hold the predicted probability
  P_OBS = 1E-10

  ; Pull the correct state probability based on the observation (DV)
  IF (DV.EQ.1) P_OBS = P_REM
  IF (DV.EQ.2) P_OBS = P_FLARE

  ; Apply the safety boundaries to the temporary variable
  IF (P_OBS.LE.0) P_OBS = 1E-10
  IF (P_OBS.GE.1) P_OBS = 0.9999999999

  ; --- ESTIMATION LOGIC ---
  Y = P_OBS 

  ; --- SIMULATION OVERWRITE ---
  ;IF (ICALL.EQ.4 .AND. EVID.EQ.0) THEN
  ;  
  ;  CALL RANDOM(2, R)
  ;  SIM_STATE = 2
  ;  
  ;  IF (R .LE. P_REM) SIM_STATE = 1
  ;  
  ;  ; Add the dummy EPS(1)*0 to satisfy NM-TRAN
  ;  Y = SIM_STATE + EPS(1)*0
  ;  
  ;ENDIF


$THETA
(0, 0.0127) ;[1] L12_BASE: Population Flaring rate
(0, 0.0656) ;[2] L21_BASE: Population Recovery rate
(0.616) ;[3] EFF_DEV: Effect of elapsed dosing intervals on flaring risk

$OMEGA
 0 FIX  ; [1] IIV rate of flaring
; uncomment the following 3 lines for simulation
;$SIGMA
;  1 FIX ; [1] Dummy EPS to satisfy NM-TRAN compiler
;$SIMULATION (12345) (54321 UNIFORM) NSUB=1000 ONLYSIM
; comment the following two lines for simulation

$ESTIMATION METHOD=1 LAPLACE LIKE PRINT=5 MAXEVALS=9999 SIGDIG=3 NOABORT
$COV

$TABLE ID TIME GA EVID DV PDV Y P_OBS P_REM P_FLARE ETA(1) PGA DOSINT GALASTD FIRSTOBS INITGA
NOAPPEND ONEHEADER NOPRINT FORMAT=tF10.5 FILE=run060.tab
