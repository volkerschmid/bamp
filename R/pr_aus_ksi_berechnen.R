pr_aus_ksi_berechnen <-
function(pr_prog, ksi_prog, noa, nop){
  
  for(i in 1:noa){
    for(j in 1:nop){
      pr_prog[j,i] <- exp(ksi_prog[j,i]) / (1+exp(ksi_prog[j,i]))
    }
  }
  
  return(pr_prog)
}
