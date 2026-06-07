pr_aus_ksi_berechnen <-
function(pr_prog, ksi_prog, noa, nop){
  pr_prog[] <- exp(ksi_prog) / (1 + exp(ksi_prog))
  return(pr_prog)
}
