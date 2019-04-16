likelihood_ratio <- function(a, b, c, d, input_type="plr", input_caller) {
  # Taken from https://stats.stackexchange.com/questions/61349/
  # Question author: Chloe T, https://stats.stackexchange.com/users/26701/
  # Response author: COOLSerdash, https://stats.stackexchange.com/users/21054/
  
  alpha <- 0.05
  
  spec <- d/(b+d)
  sens <- a/(a+c)
  
  lr.pos <- sens/(1 - spec)  
  
  if ( a != 0 & b != 0 ) {
    
    sigma2 <- (1/a) - (1/(a+c)) + (1/b) - (1/(b+d))
    lower.pos <- lr.pos * exp(-qnorm(1-(alpha/2))*sqrt(sigma2))
    upper.pos <- lr.pos * exp(qnorm(1-(alpha/2))*sqrt(sigma2)) 
    
    
  } else if ( a == 0 & b == 0 ) {
    
    lower.pos <- 0
    upper.pos <- Inf
    
  } else if ( a == 0 & b != 0 ) {
    
    a.temp <- (1/2)
    
    spec.temp <- d/(b+d)
    sens.temp <- a.temp/(a+c)
    lr.pos.temp <- sens.temp/(1 - spec.temp)  
    lower.pos <- 0
    sigma2 <- (1/a.temp) - (1/(a.temp+c)) + (1/b) - (1/(b+d))
    upper.pos <- lr.pos.temp * exp(qnorm(1-(alpha/2))*sqrt(sigma2))
    
  } else if ( a != 0 & b == 0 ) {
    
    b.temp <- (1/2)
    spec.temp <- d/(b.temp+d)
    sens.temp <- a/(a+c)
    lr.pos.temp <- sens.temp/(1 - spec.temp) 
    sigma2 <- (1/a) - (1/(a+c)) + (1/b.temp) - (1/(b.temp+d))
    lower.pos <- lr.pos.temp * exp(-qnorm(1-(alpha/2))*sqrt(sigma2))
    upper.pos <- Inf  
    
  } else if ( (a == (a+c)) & (b == (b+d)) ) {
    
    a.temp <- a - (1/2)
    b.temp <- b - (1/2)
    spec.temp <- d/(b.temp+d)
    sens.temp <- a.temp/(a+c)
    lr.pos.temp <- sens.temp/(1 - spec.temp) 
    sigma2 <- (1/a.temp) - (1/(a.temp+c)) + (1/b.temp) - (1/(b.temp+d))
    lower.pos <- lr.pos.temp * exp(-qnorm(1-(alpha/2))*sqrt(sigma2))
    upper.pos <- lr.pos.temp * exp(qnorm(1-(alpha/2))*sqrt(sigma2)) 
    
  }
  
  lr.neg <- (1 - sens)/spec
  
  if ( c != 0 & d != 0 ) {
    
    sigma2 <- (1/c) - (1/(a+c)) + (1/d) - (1/(b+d))
    lower.neg <- lr.neg * exp(-qnorm(1-(alpha/2))*sqrt(sigma2))
    upper.neg <- lr.neg * exp(qnorm(1-(alpha/2))*sqrt(sigma2)) 
    
  } else if ( c == 0 & d == 0 ) {
    
    lower.neg<- 0
    upper.neg <- Inf
    
  } else if ( c == 0 & d != 0 ) {
    
    c.temp <- (1/2)
    
    spec.temp <- d/(b+d)
    sens.temp <- a/(a+c.temp)
    lr.neg.temp <- (1 - sens.temp)/spec.temp    
    lower.neg <- 0
    sigma2 <- (1/c.temp) - (1/(a+c)) + (1/d) - (1/(b+d))
    upper.neg <- lr.neg.temp * exp(qnorm(1-(alpha/2))*sqrt(sigma2))
    
  } else if ( c != 0 & d == 0 ) {
    
    d.temp <- (1/2)
    spec.temp <- d.temp/(b+d)
    sens.temp <- a/(a+c)
    lr.neg.temp <- (1 - sens.temp)/spec.temp  
    sigma2 <- (1/c) - (1/(a+c)) + (1/d.temp) - (1/(b+d))
    lower.neg <- lr.neg.temp * exp(-qnorm(1-(alpha/2))*sqrt(sigma2))
    upper.neg <- Inf  
    
  } else if ( (c == (a+c)) & (d == (b+d)) ) {
    
    c.temp <- c - (1/2)
    d.temp <- d - (1/2)
    spec.temp <- d.temp/(b+d)
    sens.temp <- a/(a+c.temp)
    lr.neg.temp <- (1 - sens.temp)/spec.temp   
    sigma2 <- (1/c.temp) - (1/(a+c)) + (1/d.temp) - (1/(b+d))
    lower.neg <- lr.neg.temp * exp(-qnorm(1-(alpha/2))*sqrt(sigma2))
    upper.neg <- lr.neg.temp * exp(qnorm(1-(alpha/2))*sqrt(sigma2)) 
    
  }
  # Edited code from here
  if (input_type == "plr"){
    ratio <- lr.pos
    lower <- lower.pos
    upper <- upper.pos
  } else if (input_type == "nlr"){
    ratio <- lr.neg
    lower <- lower.neg
    upper <- upper.neg
  }
  
  output <- tibble(
    caller = input_caller,
    type = input_type,
    estimate = ratio,
    conf.low = lower,
    conf.high = upper,
    called = NA,
    total = NA,
    statistic = NA,
    p.value = NA,
    parameter = NA,
    method = NA,
    alternative = NA
  )
  return(output)
}
