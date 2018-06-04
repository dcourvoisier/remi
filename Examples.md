Examples
================
AU
11 mai 2018

Downloading of necessary libraries
----------------------------------

``` r
#setwd('yourworkingdirectory')
source('remi_smmr.R')
```

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

    ## Warning: package 'data.table' was built under R version 3.4.4

    ## Loading required package: Matrix

    ## Loading required package: lme4

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

``` r
library('ggplot2')
```

excitation\_function
--------------------

Random excitation pulses generation
===================================

Generation of three random pulses with amplitudes 1, 5 and 10, with duration of 2,4 and 5 seconds with measurements every second during 100 seconds, with a minimum spacing of 20 seconds between pulses.

``` r
exc <- excitation_function (amplitude = c(1,5,10), 
                            Nexc = 3, 
                            duration = c(2,4,5), 
                            deltat = 1, 
                            tmax = 100,
                            minspacing = 20)
```

    ## Warning in excitation_function(amplitude = c(1, 5, 10), Nexc = 3, duration = c(2, : Due to input parameters introduced, vector size was larger than Nf and was cut to this value.

``` r
plot(exc$t,exc$y, xlab = "Time (s)", ylab = "Excitation (unit)")
```

![](Examples_files/figure-markdown_github/unnamed-chunk-2-1.png)

calculate\_gold
---------------

Derivative calculation with uncorrelated errors
===============================================

Use of the Gold function for derivative calculation in the case of a simple quadratic function. It will be demonstrated that missing data can be present on the data and the function still manages to find the derivatives. The function chosen is:
*x*(*t*)=*t*<sup>2</sup>
 And its first and second derivatives:
$$ \\dot{x} = 2t$$
$$\\ddot{x} = 2$$

``` r
time <- c(1:500)/100
signal <- time^2
result <- calculate_gold(TimeSeries = signal, time = time, Embedding = 5)
#Puttin result data into a data table for easier plotting with ggplot
resulttable <-setDT(list(time=time,signal=signal,dtime=result$dtime,signal_rollmean=result$dsignal[,1],first_derivative=result$dsignal[,2],second_derivative=result$dsignal[,3]))
```

``` r
ggplot( data = resulttable ) +
theme_light() +
  geom_line(aes(time,signal, colour = "Signal"))+
  geom_point(aes(dtime,signal_rollmean, colour = "Signal rollmean in embedding points"))+
  geom_point(aes(dtime,first_derivative, colour = "First derivative"))+
  geom_point(aes(dtime,second_derivative, colour = "Second derivative"))+
  labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")
```

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

![](Examples_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
#Simulating missing data in input signal
time <- c(1:50)/10
time <- time[sort(sample(seq(time), 0.6*length(time),replace = F))]
signal <- time^2
result <- calculate_gold(TimeSeries = signal, time = time, Embedding = 5)
resulttable <-setDT(list(time=time,signal=signal,dtime=result$dtime,signal_rollmean=result$dsignal[,1],first_derivative=result$dsignal[,2],second_derivative=result$dsignal[,3]))
```

``` r
ggplot( data = resulttable ) +
theme_light() +
  geom_line(aes(time,signal, colour = "Signal"))+
  geom_point(aes(dtime,signal_rollmean, colour = "Signal rollmean in embedding points"))+
  geom_point(aes(dtime,first_derivative, colour = "First derivative"))+
  geom_point(aes(dtime,second_derivative, colour = "Second derivative"))+
  labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")
```

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

![](Examples_files/figure-markdown_github/unnamed-chunk-6-1.png)

remi\_generate\_order1
----------------------

Generation of the solution to the first order differential equation (convolution)
=================================================================================

``` r
exc <- excitation_function(amplitude = 10,
                            Nexc = 3, 
                            duration = 5, 
                            deltat = 0.5,
                            tmax = 100,
                            minspacing = 1)

soleq <- remi_generate_order1(dampingTime = 30,
                              inputvec = exc$y,
                              inputtim = exc$t)
excdt <- setDT(exc)
soleqdt <- setDT(soleq)
```

``` r
ggplot( ) +
theme_light() + theme(legend.position = "top") +
  geom_point(data = excdt, aes(t,y, colour = "Excitation (unit)"))+
  geom_point(data = soleqdt,aes(t,y, colour = "Convolution: solution to differential equation (unit)"))+
  labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")
```

![](Examples_files/figure-markdown_github/unnamed-chunk-8-1.png)

simulation\_generate\_order1
----------------------------

Generation of several signals for several individuals that are solution to the first order differential equation
================================================================================================================

Generating simulation data for 4 individuals, with a damping time of 10 for an excitation vector formed by 3 excitations of amplitude 1 and duration 10 s distributed randomly with a time step of 0.5 s in a total time lapse of 100 s, with a minimum spacing between pulses of 20 s and with NO NOISE. That is, the signal follows exactly the theoretical solution of the differential equation and there is no variation of the damping time, the excitation coefficient and the equilibrium value across individuals:

``` r
# Generation of signals with no noise
mydata <- simulation_generate_order1(Nindividuals = 4, 
                                    dampingTime = 10, 
                                    amplitude = 1, 
                                    Nexc = 3, 
                                    duration = 10, 
                                    deltatf = 0.5,
                                    tmax = 100,
                                    minspacing = 20,
                                    interNoise = 0, 
                                    intraNoise = 0)
```

If we add the following command we will be able to visualize the structure of mydata (the command "head" allows to visualize the first lines of the table, entering "mydata" directly will allow you to see the first and last lines).

``` r
head(mydata)
```

    ## $rawdata
    ##        ID excitation timecol Dampedsignalraw  Dampedsignal
    ##     1:  1          0    0.00   -6.389648e-15 -6.389648e-15
    ##     2:  1          0    0.01   -6.456521e-15 -6.456521e-15
    ##     3:  1          0    0.02   -6.673737e-15 -6.673737e-15
    ##     4:  1          0    0.03   -6.597116e-15 -6.597116e-15
    ##     5:  1          0    0.04   -6.030643e-15 -6.030643e-15
    ##    ---                                                    
    ## 40000:  4          0   99.96    5.441647e+00  5.441647e+00
    ## 40001:  4          0   99.97    5.436208e+00  5.436208e+00
    ## 40002:  4          0   99.98    5.430774e+00  5.430774e+00
    ## 40003:  4          0   99.99    5.425346e+00  5.425346e+00
    ## 40004:  4          0  100.00    5.419923e+00  5.419923e+00
    ## 
    ## $data
    ##      ID excitation timecol Dampedsignalraw  Dampedsignal
    ##   1:  1          0     0.0   -6.389648e-15 -6.389648e-15
    ##   2:  1          0     0.5   -1.858133e-15 -1.858133e-15
    ##   3:  1          0     1.0    2.251170e-15  2.251170e-15
    ##   4:  1          0     1.5    2.501739e-15  2.501739e-15
    ##   5:  1          0     2.0    2.316254e-16  2.316254e-16
    ##  ---                                                    
    ## 800:  4          1    98.0    6.569759e+00  6.569759e+00
    ## 801:  4          0    98.5    6.297053e+00  6.297053e+00
    ## 802:  4          0    99.0    5.989942e+00  5.989942e+00
    ## 803:  4          0    99.5    5.697809e+00  5.697809e+00
    ## 804:  4          0   100.0    5.419923e+00  5.419923e+00

Where: ID is the identifier of the individual excitation is the excitation signal Dampedsignalraw is the signal without noise Dampedsignal is the signal with noise timecol is the time column generated.

Plotting data:

``` r
ggplot( data = mydata$data ) +
  geom_point(aes(timecol,Dampedsignalraw, colour = "Signal-no noise"))+
  geom_point(aes(timecol,Dampedsignal, colour = "Signal with 0% intra-noise"))+
  geom_line(aes(timecol,excitation,colour = "Excitation"))+
  facet_wrap(~ID)+
  labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")
```

![](Examples_files/figure-markdown_github/unnamed-chunk-11-1.png)

Using the same function, this time adding a 20% intra-individual noise and a 40% inter-individual noise:

``` r
# Generation of signals with intra and inter-noise
mydata <- simulation_generate_order1(Nindividuals = 4, 
                                    dampingTime = 10, 
                                    amplitude = 1, 
                                    Nexc = 3, 
                                    duration = 10, 
                                    deltatf = 0.5,
                                    tmax = 100,
                                    minspacing = 20,
                                    interNoise = 0.4, 
                                    intraNoise = 0.2)
```

    ## Warning in excitation_function(amplitude, Nexc, duration, deltat, tmax, : Due to input parameters introduced, vector size was larger than Nf and was cut to this value.

    ## Warning in excitation_function(amplitude, Nexc, duration, deltat, tmax, : Due to input parameters introduced, vector size was larger than Nf and was cut to this value.

Plotting data:

``` r
ggplot( data = mydata$data ) +
  geom_point(aes(timecol,Dampedsignalraw, colour = "Signal-no noise"))+
  geom_point(aes(timecol,Dampedsignal, colour = "Signal with 20% intra-noise"))+
  geom_line(aes(timecol,excitation,colour = "Excitation"))+
  facet_wrap(~ID)+
  labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")
```

![](Examples_files/figure-markdown_github/unnamed-chunk-13-1.png)

remi\_analyse\_order1
---------------------

Study of input signals to indicate how well it fits a first order differential equation through multilevel regression
=====================================================================================================================

Next, the signals with noise presented above will be analyzed in order to verify that the damping coefficient was the one introduced in the simulation function and that the estimated signal generated matches the original one:

``` r
result <- remi_analyse_order1(UserData = mydata$data,
                                ID = "ID",
                                Input="excitation",
                                Time="timecol",
                                signalcolumn = "Dampedsignal",
                                Embedding = 5)
```

    ## summary from lme4 is returned
    ## some computational error has occurred in lmerTest

Now let's take a look to the different parts of the result. As it was mentioned before, the first table contains the original data with some columns added. These columns contain intermediate variables used for the preparation of the regression:

``` r
head(result$data)
```

    ##    ID excitation timecol Dampedsignalraw Dampedsignal
    ## 1:  1          0     0.0   -7.539647e-15    1.9893913
    ## 2:  1          0     0.5   -3.911356e-15   -1.0099792
    ## 3:  1          0     1.0   -9.596821e-16   -3.3340040
    ## 4:  1          0     1.5    2.400093e-15    2.1662179
    ## 5:  1          0     2.0   -4.008284e-15   -1.6519624
    ## 6:  1          0     2.5   -1.982024e-15   -0.2960264
    ##    Dampedsignal_rollmean Dampedsignal_derivate1 timecol_derivate
    ## 1:            -0.3680673             -0.8213020              1.0
    ## 2:            -0.8251508              0.6219894              1.5
    ## 3:            -0.5341579              1.0191469              2.0
    ## 4:            -0.2216410             -1.1556655              2.5
    ## 5:            -1.0696692             -0.4638628              3.0
    ## 6:            -0.7222585             -0.3513348              3.5
    ##    excitation_rollmean Dampedsignal_estimated
    ## 1:                   0            -0.03047207
    ## 2:                   0            -0.03047207
    ## 3:                   0            -0.03047207
    ## 4:                   0            -0.03047207
    ## 5:                   0            -0.03047207
    ## 6:                   0            -0.03047207

Where: Dampedsignal\_rollmean contains the roll mean (moving average) values of the input signal in embedding points. As it can be seen, the first line contains an NA because the convolution takes the points to the left and thus the first roll means can't be calculated as there are no points to the left of these.

Dampedsignal\_derivate1 contains the first derivate of Dampedsignal, calculated by using the calculate\_gold function. The first line contains an NA for the same reason as the previous column.

timecol\_derivate contains the values of time in which the derivative has been evaluated.

excitation\_rolled contains the roll mean of the excitation signal in embedding points.

Dampedsignal\_estimated contains the values of the estimated signal generated by using the remi\_generate\_order1 function and using the coefficients calculated for each individual (see next table).

``` r
result$resultID
```

    ##    ID Dampedsignal_dampingTime Dampedsignal_eqvalue excitation_exccoeff
    ## 1:  1                 8.998980          -0.03047207            8.914529
    ## 2:  2                 8.999026          -0.03047222            8.914581
    ## 3:  3                13.426478          -0.04546432           13.794718
    ## 4:  4                11.976019          -0.04055282           12.195944

Where for each individual we have:

Excitation\_exccoeff which is the coefficient of the excitation term.

Dampedsignal\_dampingTime which is the inverse of the damping coefficient.

Dampedsignal\_eqvalue which is the equilibrium value.

(for more details, visit the wiki pages)

Finally, the third table contains the average of these coefficients for all the individuals:

``` r
result$resultmean
```

    ##     ID Dampedsignal_dampingTime Dampedsignal_eqvalue excitation_exccoeff
    ## 1: All                 10.52004          -0.03562263            10.59111
    ##    Dampedsignal_fitmsg
    ## 1:                <NA>

And, if we want to see in detail the results of the regression, the following command applies:

``` r
result$regression
```

    ## [[1]]
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## Dampedsignal_derivate1 ~ Dampedsignal_rollmean + excitation_rollmean +  
    ##     (1 + excitation_rollmean + Dampedsignal_rollmean | ID)
    ##    Data: Data
    ## Control: lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap")
    ## 
    ## REML criterion at convergence: 1963.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5718 -0.6951 -0.0355  0.6749  3.1769 
    ## 
    ## Random effects:
    ##  Groups   Name                  Variance  Std.Dev. Corr     
    ##  ID       (Intercept)           0.0000000 0.00000           
    ##           excitation_rollmean   0.0005386 0.02321   NaN     
    ##           Dampedsignal_rollmean 0.0005337 0.02310   NaN 1.00
    ##  Residual                       0.6936782 0.83287           
    ## Number of obs: 788, groups:  ID, 4
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)           -0.003386   0.047045  -0.072
    ## Dampedsignal_rollmean -0.095057   0.019432  -4.892
    ## excitation_rollmean    1.006755   0.078646  12.801
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Dmpds_
    ## Dmpdsgnl_rl -0.508       
    ## excttn_rllm -0.090 -0.288
    ## 
    ## [[2]]
    ## $ID
    ##   (Intercept) excitation_rollmean Dampedsignal_rollmean
    ## 1           0         -0.01614003           -0.01606702
    ## 2           0         -0.01613929           -0.01606645
    ## 3           0          0.02067095            0.02057699
    ## 4           0          0.01160837            0.01155648

Where we have a summary of the random and fixed effects and the residuals calculated by the function lmer. Apart from these indicators, if we graphically wish to verify how the estimated signal fits the initial signal, we can call ggplot once again:

``` r
ggplot( data = result$data ) +
  geom_point(aes(timecol,Dampedsignal_estimated, colour = "Estimated signal"))+
  geom_line(aes(timecol,Dampedsignalraw, colour = "Signal-no noise"))+
  geom_point(aes(timecol,Dampedsignal, colour = "Signal-20% intra-noise"))+
  geom_line(aes(timecol,excitation,colour = "Excitation"))+
  facet_wrap(~ID)+
  labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")
```

![](Examples_files/figure-markdown_github/unnamed-chunk-19-1.png)
