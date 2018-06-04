#Necessary libraries
library(zoo) # for the calculation of rollmeans
library(data.table) # for the data processing syntax
library(lmerTest) # for the calculation of p values for regressions
library(lme4) # for the calculation of regressions


# calculate_gold ----------------------------------------------------------


# GOLD function for derivative calculation
calculate_gold <-  function(TimeSeries,
                            time,
                            Embedding=2)
{
  
  #Error management
  if (length(TimeSeries)!=length(time)) {
    stop("TimeSeries and time vectors should have the same length.\n")
  }
  if(length(TimeSeries)<=Embedding) {
    stop("TimeSeries and time vectors should have a length greater than Embedding. Not enough points present to do the calculations.Change TimeSeries and time vectors or reduce the Embedding value.\n")
  }
  if(Embedding > 2){
    #tembed (time embedded) is a matrix, containing in each line the groups of "n" points with n being the embedding value from which to calculate the roll means
    #For instance if time=1:10 and Embedding=3, the matrix tembed will have: row1: 1 2 3; row2: 2 3 4...row8:8 9 10).
    tembed <- embed(time,Embedding) #The "stats" library "embed" function provides a matrix in which the groups of values are inversed
    #by column and thus the next operation is to format these values so that the matrix contain the groups by line and from left to right.
    tembed<-tembed[,ncol(tembed):1]
    
    #Creation of the D matrix for estimation of derivatives up to second derivative.
    #According to the paper mentioned, equations (11) and (12) this is "...the diagonal matrix with scaling constants that convert the
    #polynomial estimates to derivative estimates". 
    D <- cbind(c(1,0,0),c(0,1,0),c(0,0,0.5)) #Cbind does a column binding of the vectors
    
    
    #Creation of the empty derivative matrix
    #Lines: It will have as many lines as groups of time points the embed funcion gives (thus, the number of lines of the tembed matrix)
    #Columns: It will have 3 columns because:
    #   First column will contain the signal mean value in the time points considered in "Embedding"
    #   Second column will contain the signal derivative in those same time points 
    #   Third column will contain the signal second derivative in those same time points
    derivative <- matrix(NA,nrow = nrow(tembed),ncol = 3) 
    
    
    #Xembed is a matrix containing the values of the signal in the time points considered in "Embedding" (no mean calculated yet)
    #As for tembed, it contains in each line the groups of "n" points of the function with n being the embedding value to PREPARE for the calculation of the roll mean
    Xembed = embed(TimeSeries,Embedding)
    Xembed<-Xembed[,ncol(Xembed):1]
    
    for(k in 1:nrow(tembed)) # Loop repeated in each tembed line (each group of time values)
      #To take into account variable delta t.
    {
      
      t <- tembed[k,]-tembed[k,(Embedding+1)/2] #Time vector, containing deltat, centered in 0
      E <- matrix(NA,nrow = 3, ncol = length(t)) #Construction of the empty Theta matrix, equation (10) 
      #But will be filled from the polynomials calculated in equation (9). 
      #It has 3 rows as we are calculating up to the second derivative.
      
      #The two for loops resolve the paper's equation (9) (polynomial's coefficients)
      for(i in 1:length(t))
      {
        E[1,i] = 1
        E[2,i] = t[i]
        
      }
      for(i in 1:length(t))
      {
        E[3,i] = t[i]^2-sum(t^2)/(sum(E[1,]))-t[i]*sum(t^3)/sum(t^2)
      }
      
      # And with E (Theta) and D it is possible to calculate the orthogonal matrix W, equation (14)
      L <- D%*%E
      W <- t(L)%*%solve(L%*%t(L))
      
      #And with this matrix, solve the differential equation as Y=X*W
      derivative[k,] <- Xembed[k,] %*% W # Derivative calculated in the time row k
    }
    derivative <- rbind(derivative,matrix(data = NA, ncol = 3,nrow = Embedding-1)) # Addition of NA so that the derivative rows are the same as those of Timeserie
    
  } else if(Embedding == 2){
    warning("Only first derivative can be calculated with an Embedding of 2.\n")
    derivative <- cbind(rollmean(TimeSeries,Embedding), diff(TimeSeries)/diff(time)) 
    #Appends the two columns: 1. mean time values and 2. The span calculated as the difference of two signal values (going forward) divided by the time interval
    derivative <- rbind(derivative,matrix(data = NA, ncol = 2,nrow = Embedding-1))
    # Addition of NA so that the derivative rows are the same as those of Timeseries
  } else{
    stop("Embedding should be at least 2 for the calculation of a first derivative and at least 3 for the calculation of a second derivative.\n")
  }
  time_derivative <-c(rollmean(time,Embedding),rep(NA,Embedding-1)) # Addition of NA so that the time rows are the same as those of Timeseries
  
  returnobject <- list("dtime" = time_derivative, 
                       "dsignal" = derivative,
                       "Embedding" = Embedding)
  return(returnobject)
}

# excitation_function -----------------------------------------------------


# Excitation signal generation.
excitation_function = function(amplitude = 1,
                               Nexc = 1,
                               duration = 2,
                               deltat = 0.1,
                               tmax = 10,
                               minspacing = 1)
{
  #Error management
  if (any(duration<=0)|any(amplitude==0)|Nexc<=0) {
    stop("Invalid input parameters. At least one excitation must be defined. Duration and Nexc must be greater than 0.Amplitude must be different from 0.\n")
  }
  if (Nexc<length(duration)|Nexc<length(amplitude)) {
    warning("The number of excitations Nexc is smaller than the number of elements in amplitude and duration. Only the first elements of these vectors were considered.\n")
  }
  
  if (Nexc>length(duration) && length(duration)>1) {
    duration<-rep(duration,ceiling(Nexc/length(duration)))
    duration<-duration[1:Nexc]
    warning("The number of excitations Nexc was higher than the durations defined. The values from the vector were repeated.\n")
  }
  if (Nexc>length(amplitude) && length(amplitude)>1) {
    amplitude<-rep(amplitude,ceiling(Nexc/length(amplitude)));
    amplitude<-amplitude[1:Nexc]
    warning("The number of excitations Nexc was higher than the amplitudes defined. The values from the vector were repeated.\n")
  }

  if(tmax<(sum(duration)+minspacing)*Nexc){stop("Non valid parameters. tmax should be greater than (duration+minspacing)*Nexc.\n")}
  
  if (length(duration)==1) {dur<-duration*Nexc} #If it's a scalar, it assigns the value directly. If it is not, it takes the value from the vector
  #at the position of the correponding excitation
  else {dur<-sum(duration)}
  
  trest<-tmax-dur-minspacing*(Nexc-1) #at the beginning, calculates the time left (to calculate random time from it) is tmax-total pulse duration - total minspacing duration
  cumt<-0 #cumulated time initialized to 0
  
  for (i in 1:Nexc)
  {
    #Generation of a single unitary pulse (single pulse with amplitude=1)
    if (length(duration)==1) {dur<-duration;durleft<-dur*(Nexc-i)} #If it's a scalar, it assigns the value directly. If it is not, it takes the value from the vector
    #at the position of the correponding excitation
    else {dur<-duration[i];durleft<-sum(duration[(i+1):Nexc])}
    
    tal<-tmax #Random time is initialized so that it enters the while loop the first time
    
    while(tal>trest){ #While the random time coming from the sample is greater than the time left before tmax
      if(i==1){tal<-sample(0:trest,1,replace=T)}
      else{
        if(trest==minspacing){tal<-minspacing}
        else{tal<-sample(minspacing:trest,1,replace=T)}
      }
    }
    Nf<-tmax/deltat+1
    sp<-rep(c(0,1),c(tal/deltat+1,dur/deltat+1)) #simple unitary pulse
    cumt<-cumt+tal+dur #cumulated time
    
    #Generation of amplified simple pulse
    if (length(amplitude)==1) {amp<-amplitude} #Same as for duration.
    else {amp<-amplitude[i]}
    sp<-sp*amp
    
    #Append calculated single pulse to result vector
    if(i==1){E<-sp}
    else{E<-append(E,sp)}
    
    #Calculating the time left so that a new sample of random time can be taken
    trest<-tmax-cumt-durleft-minspacing*(Nexc-i-1)
    
  }
  
  #Verify E vector length
  if (length(E)<Nf) {
    E<-append(E,rep(0,Nf-length(E))) #If the length is smaller than Nf, fill with 0
  } 
  else{
    E<-E[1:Nf]
    warning("Due to input parameters introduced, vector size was larger than Nf and was cut to this value.\n")
  }
  #Generation of time vector
  tim<-seq(0,tmax,deltat)
  
  data<-list(y=E,t=tim)
  return(data)    
  
}


# remi_analyse_order1 -----------------------------------------------------


# remi first order analysis function
remi_analyse_order1 <- function(UserData, 
                                ID = NULL, 
                                Input = NULL, 
                                Time = NULL, 
                                signalcolumn,   
                                Embedding = 2){
  
  Data<-copy(UserData) #Keeps the original of Data so that it can rename columns freely
  #Taken from: https://stackoverflow.com/questions/15913417/why-does-data-table-update-namesdt-by-reference-even-if-i-assign-to-another-v
  
  Data <- setDT(Data) # Convert input data to data.table. 
  
  noinput<-FALSE #Flag that will allow to differentiate if there is an excitation term or not when doing the regression
  #Error management
  #Verifying column names repeated in data table.
  if(any(duplicated(colnames(Data)))){stop("Input datatable contains duplicated column names. Please correct these in order to launch the analysis function.\n")}
  
  #If ID,Input,Time,signalcolumn are not strings containing the name of columns in "Data", stop the function
  if(is.null(ID)|(!is.null(ID) && !is.character(ID))){stop("ID should be a string containing the name of the column in Data that contains the individual identifier.\n")}
  if(!is.null(Input) && !is.character(Input)){stop("Input should be a string containing the name of the column in Data that contains the excitation.\n")}
  if(!is.null(Time) && !is.character(Time)){stop("Time should be a string containing the name of the column in Data that contains the time.\n")}
  if(!is.character(signalcolumn)){stop("signalcolumn should be a string containing the name of the column in Data that contains the signal of the individual.\n")}
  
  if(is.null(Input)){Data[,Inputcol := 0]
    Input = "Inputcol"
    noinput<-TRUE #This flag will be needed later as if there is no input, coefficients for the excitation term will not be calculated in the regression 
    warning("No excitation signal introduced as input. Input is set to 0.\n")}
  
  if(is.null(Time)){Data[,timecol := c(1:.N),by = ID]
    Time <- "timecol" # if no time set it to a 1 sec step vector
    warning("No time vector introduced as input. A 1 unit increment time vector is generated.\n")}

    #Rename column containing ID to "ID" because further functions use this identifier to work. Also display warning message saying that column has been renamed.
  if(ID!="ID"){
    #Setting some of the column names (necessary for data treatment)
    setnames(Data, old=ID, new="ID")}
  
  if(Embedding==2){warning("Only first derivative can be calculated with an Embedding of 2.\n")}
  
  #Suppress NA elements if there is an NA in Time or in signalcolumn
  #This suppresses the entire row of the data table
  Data <- Data[!is.na(get(Time)) & !is.na(get(signalcolumn))] 
  
  #If in the data left, there are some NA values in the excitation vector, set them to zero.
  #This is to avoid loosing data from signal and time vectors, as sometimes when people fill in the signal data they put NA to mean no excitation, thus 0.
  Data[is.na(get(Input))]<-0
  
  #Calculation of the signal rollmean and first derivative of the signal column
  #Paste is only used to generate new column names based on the orignal ones (concatenate strings)
  options(warn = -1) #Turns off warnings to avoid the calculate_gold function activating warnings every time it is called
  
  Data[,c(paste0(signalcolumn,"_rollmean")) := calculate_gold(get(signalcolumn),get(Time),Embedding)$dsignal[,1], by = ID]
  Data[,c(paste0(signalcolumn,"_derivate1")) := calculate_gold(get(signalcolumn),get(Time),Embedding)$dsignal[,2], by = ID]
  Data[,c(paste0(Time,"_derivate")) := calculate_gold(get(signalcolumn),get(Time),Embedding)$dtime, by = ID]
  
  options(warn = 0) #Turns warnings back on
  
  #Calculation of the roll mean of the excitation column
  Data[ ,(paste0(Input,"_rollmean")) := c(rollmean(get(Input),Embedding),rep(NA,Embedding-1)), by=ID] 
 
  #Result data frames
  #The first table contains the input data and the calculations made above
  
  #The second table contains the mean results for excitation coefficient, damping time and offset for each individual
  ResultID <- setDT(list(ID = unique(Data$ID)))
  
  #The third table contains the mean values for fixed part of coefficients (divided by damping time in the case of the excitation and offset) for all the individuals (single line) 
  
  # First order derivative equation mixed regression 
  setkey(ResultID,ID) #sorts the data table by ID
    
    if(noinput){ # if there is no excitation signal
      model <- tryCatch({lmer(paste0(signalcolumn,"_derivate1 ~ ",signalcolumn,"_rollmean + (1 + ",signalcolumn,"_rollmean |ID)"),
                               data=Data, REML=TRUE,
                               control = lmerControl(calc.derivs = FALSE,optimizer = "nloptwrap"))}, error=function(e) e)
    }else{ # if there is an excitation signal
      model <- tryCatch({ lmer(paste0(signalcolumn,"_derivate1 ~ ",signalcolumn,"_rollmean + ",Input,"_rollmean + (1 + ",Input,"_rollmean + ",signalcolumn,"_rollmean |ID)"),
                               data=Data, REML=TRUE,
                               control = lmerControl(calc.derivs = FALSE,optimizer = "nloptwrap"))}, error=function(e) e)}
    if (!inherits(model,"error")){ # if the regression worked
      summary <- summary(model) # Summary of the regression
      random <- ranef(model) # Variation of the estimate coefficient over the individuals
      regression <- list(summary, random) # list to output both results: summary, and the table from ranef
      
      #Generate mean results (third output table) with convergence criterions
      
      # Extract the damping coefficient from the regression summary
      Resultmean<-setDT(list(ID="All"))
     
      # calculate the damping time for all signal columns
      Resultmean[, c(paste0(signalcolumn,"_dampingTime")) :=  
                   -1L/summary$coefficients[paste0(signalcolumn,"_rollmean"),"Estimate"]] # calculate the damping time: -1/damping_coeff
      
      # Extract the intercept coeff (equilibrium value)
      Resultmean[,c(paste0(signalcolumn,"_eqvalue")) := summary$coefficients["(Intercept)","Estimate"]*Resultmean[, get(paste0(signalcolumn,"_dampingTime"))]]
     
      #Generate the results for each individual (second output table)
      #Calculate the damping time from the damping coeff
      
      # Damping time in ResultID ------------------------------------------------
      
      
      ResultID[, c(paste0(signalcolumn,"_dampingTime")) := -1L/(summary$coefficients[paste0(signalcolumn,"_rollmean"),"Estimate"] + random$ID[.GRP,paste0(signalcolumn,"_rollmean")]), by = ID] 
      
      
      # Extract the intercept (equilbrium value) calculated for each individual (present in random, regression table)
      # Offset in ResultID ------------------------------------------------------
      ResultID[,c(paste0(signalcolumn,"_eqvalue")) := (summary$coefficients["(Intercept)","Estimate"]+ random$ID[.GRP,"(Intercept)"])*ResultID[.GRP, get(paste0(signalcolumn,"_dampingTime"))], by=ID] 
      
      # Generation of the fitted signal for all ID using remi generate 
      if(!noinput){# If there is an excitation term
        # Extract the excitation coeff for the excitation
          Resultmean[ ,c(paste0(Input,"_exccoeff")) := summary$coefficients[paste0(Input,"_rollmean"),"Estimate"]*Resultmean[, get(paste0(signalcolumn,"_dampingTime"))]]

          #And for each individual: the mean coeff (sumary$coeff) + the variation per Individual (in random)
          # Excitation coefficient in ResultID --------------------------------------
          ResultID[,c(paste0(Input,"_exccoeff")) :=
                     (summary$coefficients[paste0(Input,"_rollmean"),"Estimate"] + random$ID[.GRP,paste0(Input,"_rollmean")])*ResultID[.GRP, get(paste0(signalcolumn,"_dampingTime"))], by = ID]}
      
          Data[,c(paste0(signalcolumn,"_estimated")) :=
             if(!is.na(ResultID[.GRP, get(paste0(signalcolumn,"_dampingTime"))]) && ResultID[.GRP, get(paste0(signalcolumn,"_dampingTime"))]>0 ){
               # if there is a damping time that has been calculated and if it is greater than 0 (decreasing exponential)
               #and if there is an excitation coefficient (excitation term was found in the inputs)
               if(!noinput){ #There is an excitation signal as input
                 remi_generate_order1(ResultID[.GRP, get(paste0(signalcolumn,"_dampingTime"))],get(Input)*(summary$coefficients[paste0(Input,"_rollmean"),"Estimate"] + random$ID[.GRP,paste0(Input,"_rollmean")]),get(Time))$y+
                 ResultID[.GRP, get(paste0(signalcolumn,"_eqvalue"))]}
               else{ #There is no excitation signal as input
               "There is no excitation provided, thus signal wasn't generated."}}
             else{NaN}, by = ID]}
  
    else { # if the regression didn't work, set to NA all coeffs
      
      ResultID[,c(paste0(signalcolumn,c("_dampingTime","_exccoeff","_eqvalue"))) := NA]
      Resultmean[,c(paste0(signalcolumn,c("_dampingTime","_exccoeff","_eqvalue"))) := NA]
      regression <- model}

    # Write any error message from the regression
    Resultmean[,c(paste0(signalcolumn,"_fitmsg")) := summary$fitMsgs]
    # Output the results for the function
    return(list(data = Data, resultID = ResultID, resultmean = Resultmean, regression = regression))}


# remi_generate_order1 ----------------------------------------------------
#Generation of the solution to the first order differential equation
remi_generate_order1 = function(dampingTime,
                                inputvec,
                                inputtim)
{
  #Error management
  #If inputvec is a scalar, the function warns the user that it should be a vector containing the values of the excitation signal
  if (length(inputvec)<=1|length(inputvec)!=length(inputtim)) {
    warning("Both the excitation (inputvec) and is time values (inputtim) should be vectors and have the same length.\n")
  }
  #if inputvec is a character or a matrix, the function stops
  if (is.matrix(inputvec) | is.character(inputvec)) stop("Inputvec should be a vector.")
  
  green <- exp(-1/dampingTime*(inputtim-inputtim[1])) #Calculation of green function values.Inputtim vector is forced to start in 0
  
  #as convolution doesn't take time into account, only number of points, a coefficient is needed to reduce amplitude
  convol <- (max(inputtim)-min(inputtim))/(length(inputtim)-1)*rev(convolve(rev(green),inputvec,type = "open")) 
  convol <- convol[1:length(inputvec)]
  return(list(y=convol,t=inputtim))              #Returns the result of the convolution in the points corresponding to the excitation time vector.
}


# simulation_generate_order1 ----------------------------------------------

# Simulation of various individual signals with intra and inter noise
simulation_generate_order1 = function(Nindividuals = 1, 
                                      dampingTime, 
                                      amplitude = 1, 
                                      Nexc = 1,
                                      duration = 10, 
                                      deltatf=1,
                                      tmax=10,
                                      minspacing = 10, 
                                      interNoise = 0,
                                      intraNoise = 0)
{
  #Assuming that dampingTime/deltat>=10.
  if (dampingTime<10) deltat=dampingTime/10
  else deltat=0.01 #Internal parameter to generate pseudo-continuous function
 
  if((deltatf %% deltat) != 0){stop("Invalid deltatf, please modify.\n")}
 
  Npoints<-tmax/deltat+1
  # Generate simulation data for a given excitation and damping time
  # Generates ID column (ID being the individual number)
  # Creates a data table with ID being the first column containing as many lines per individual as time points marked in Npoints
  
  Data <- setDT(list(ID = unlist(lapply(c(1:Nindividuals),function(x){rep(x,Npoints)})))) 
  
  #Add to that data table an "excitation" column, containing the values of the excitation signal. 
  #Creates a new excitation signal for each individual
  Data[,excitation := excitation_function(amplitude,Nexc,duration,deltat,tmax,minspacing)$y,by = ID]
  Data[,timecol := excitation_function(amplitude,Nexc,duration,deltat,tmax,minspacing)$t,by = ID] 
  
  #Creates a damping time vector by taking dampingTime input value and adding the interNoise in a normal distribution
  #Contains as many elements as individuals
  dampingTimevec = dampingTime + rnorm(Nindividuals, mean = 0, sd = interNoise*dampingTime)
  
  #If any value of the damping time vector is negative, the original value is used instead at that position of the vector
  if(any(dampingTimevec<0)){
    dampingTimevec[dampingTimevec < 0] <- dampingTime  
    warning("Some values for dampingTime where negative when adding interNoise. The original dampingTime value has been used instead for those cases")
  }
  
  #Creates the signals for each individual taking the damping time for that individual from dampingTimevec
  #Dampedsignalraw is the signal WITHOUT NOISE
  Data[,Dampedsignalraw := remi_generate_order1 (dampingTimevec[.GRP],excitation,timecol)$y, by = ID ] 
  
  #Creates the signal for each individual with intra noise
  #In order to avoid adding an increased intra noise that will "grow" with each new pulse of the excitation signal, first
  #a "amplitudenorm" function is calculated in which the excitation has only one pulse that starts at the begining of the time sequence
  #The intranoise is added regarding the maximum value of that function
  #As the excitation function can receive an amplitude vector and a duration vector, in order to normalize, it will
  #be necessary to pick the maximum of the amplitude and the excitation
  if(length(amplitude)>1){amp<-max(amplitude)}
  else {amp<-amplitude}
  if(length(duration)>1){dur<-max(duration)/deltat+1} #As done in the excitation function
  else{dur<-duration/deltat+1}
  Data[,amplitudenorm := remi_generate_order1(dampingTimevec[.GRP],rep(c(amp,0),c(dur,(length(excitation)-dur))),timecol)$y, by = ID ]
  Data[,Dampedsignal := Dampedsignalraw + rnorm(.N, mean = 0, sd = intraNoise*max(abs(amplitudenorm))), by = ID ] 
  
  #Returning rawdata
  RawData<-copy(Data)
  
  #Selecting equally spaced data with deltatf
  #definition of time step for equally spaced values. 
  Data<-Data[timecol %in% seq(0,tmax,deltatf), .SD,by = ID] #Keeping only values from data table that are in the time intervals chosen
  return(list(rawdata=RawData[,!"amplitudenorm",with = FALSE],data=Data[,!c("amplitudenorm"),with = FALSE]))
}
