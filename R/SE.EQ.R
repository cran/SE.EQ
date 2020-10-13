SE.EQ <-
function( X , Y , eq_margin , alpha=0.05 , print.results=TRUE){
  
  if (print.results) {
    cat("\n")
    cat("*******************************************************************************" , "\n\n")
    cat("********************* \t The SE-test for equivalence \t***********************" , "\n")
    cat("*************** \t ( see Hoffelder et al. (2015) \t***********************"      , "\n\n")
    cat("*******************************************************************************" , "\n\n")
  }
  
  m_X         <- nrow(X)                                                # sample size first sample
  n_Y         <- nrow(Y)                                                # sample size second sample
  N           <- m_X + n_Y                                              # total sample size pooled sample
  d           <- ncol(X)                                                # number of variables / dimension
  Mean_X      <- colMeans(X)                                            # Mean of first sample
  Mean_Y      <- colMeans(Y)                                            # Mean of second sample 
  Mean_diff   <- Mean_X - Mean_Y                                        # Difference of both sample mean vectors
  S_X         <- var(X)                                                 # empirical covariance matrix of first sample
  S_Y         <- var(Y)                                                 # empirical covariance matrix of second sample
  vars_X      <- diag(S_X)                                              # empirical variances of first sample
  vars_Y      <- diag(S_Y)                                              # empirical variances of second sample
  pi1hat      <- 2 * S_X * S_X                                          # Matrix pi1hat according to Hoffelder et al. (2015)
  pi2hat      <- 2 * S_Y * S_Y                                          # Matrix pi2hat according to Hoffelder et al. (2015)
  
  Q_SE        <- sum( (Mean_diff * Mean_diff) / (0.5 * (vars_X + vars_Y) ) )   # estimated distance of SE-test

  # calc of asymptotic variance
  MATRIX_HELP_1 <-  S_X/m_X + S_Y/n_Y
  MATRIX_HELP_2 <-  pi1hat/m_X + pi2hat/n_Y
  delt_hat    <-  Mean_diff / (vars_X + vars_Y)
  eps_hat     <-  delt_hat * delt_hat
  asympt_var  <-  16*t(delt_hat) %*% MATRIX_HELP_1 %*% delt_hat + 4*t(eps_hat) %*% MATRIX_HELP_2 %*% eps_hat          # estimated asymptotic variance
 # cat("vars_X",vars_X,"\n")
#  cat("vars_Y",vars_Y,"\n")
#  cat("pi1hat",pi1hat,"\n")
#  cat("pi2hat",pi2hat,"\n")
#  cat("MATRIX_HELP_1",MATRIX_HELP_1,"\n")
# cat("MATRIX_HELP_2",MATRIX_HELP_2,"\n")
# cat("delt_hat",delt_hat,"\n")
# cat("eps_hat",eps_hat,"\n")
  numerator   <- Q_SE -  eq_margin                                    # numerator of SE test statistic
  denominator <- sqrt(asympt_var)                                     # denominator of SE test statistic

  quantile		<- qnorm( p=alpha , mean = 0, sd = 1 )                    # quantile of standard normal distribution
  teststat		<- numerator / denominator 	                              # teststatistic
  
  if (print.results) {
    cat(" Summary statistics:"                                                           , "\n\n",    
        "Sample size of first sample X:"                   ,"\t" ,m_X                    ,   "\n",    
        "Sample size of second sample Y:"                  ,"\t" ,n_Y                    ,   "\n",
        "Sample size of pooled sample (X,Y):"              ,"\t" ,N                      ,   "\n",
        "Dimension (Number of variables):"                 ,"\t" ,d                      ,   "\n",
        "Mean of first sample:"                        ,"\t\t\t" ,Mean_X                 ,   "\n",
        "Mean of second sample:"                         ,"\t\t" ,Mean_Y                 ,   "\n",
        "Difference of both means:"                      ,"\t\t" ,Mean_diff              ,   "\n",
        "Empirical covariance matrix S_X:"                 ,"\t" ,round(S_X,digits=2)    ,   "\n",
        "Empirical covariance matrix S_Y:"                 ,"\t" ,round(S_Y,digits=2)    ,   "\n",
        "Estimated SE distance:"                        ,"\t\t"  ,Q_SE                   ,   "\n",      
        "Equivalence margin:"                          ,"\t\t\t" ,eq_margin              ,   "\n",      
        "asymptotic variance:"                         ,"\t\t\t" ,asympt_var             ,   "\n",   
        "Significance level:"                          ,"\t\t\t" ,alpha                  ,   "\n"
    )
  }
  
  
  if 	( teststat < quantile ) 	{
    erg_text	<-	'Equivalence comparison successful'	                  # Test result text 
    erg       <- 1                                                      # Test result  
  }
  else 	   {
    erg_text	<-	'Equivalence comparison not successful'		            # Test result text 
    erg       <- 0                                                      # Test result  
  }
  
  p_value			<- pnorm( q=teststat , mean = 0, sd = 1 )                 # p-value of SE-test for equivalence  
  
  if (print.results) {
    cat(" Teststatistic:"                                        ,"\t\t\t" ,teststat ,   "\n",    
        "Quantile of stand. normal dist.:"                       ,"\t"     ,quantile ,   "\n\n",
        "Decision in favor (1) or against (0) equivalence: "               ,erg      ,   "\n\n",
        "Test result:"                                                               ,   "\n\n",
        "\t p-value of the SE-test for equivalence: p ="                   ,p_value  ,   "\n\n",
        "\t\t"                                                             ,erg_text ,   "\n\n",
        "******************************************************************************"  ,   "\n\n"
    )
  }
  result.summary <- data.frame(p.value = p_value , testresult.num = erg , testresult.text = erg_text )
  result.summary	          # Return of Test results  
}
