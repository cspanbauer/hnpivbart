## combine_bart <- function(bart.list){
##     chains <- length(bart.list)
##     out <- bart.list[[1]]
##     out$chain <- rep(1,nrow(out$dmu1))
##     for(chn in 2:chains){
##         out$dnpart <- c(out$dnpart,bart.list[[chn]]$dnpart)
##         out$dalpha <- c(out$dalpha,bart.list[[chn]]$dalpha)
##         out$dmu1 <- rbind(out$dmu1,bart.list[[chn]]$dmu1)
##         out$dsigma1 <- rbind(out$dsigma1,bart.list[[chn]]$dsigma1)
##         out$dmu2 <- rbind(out$dmu2,bart.list[[chn]]$dmu2)
##         out$dsigma2 <- rbind(out$dsigma2,bart.list[[chn]]$dsigma2)
##         out$dgamma <- rbind(out$dgamma,bart.list[[chn]]$dgamma)
##         out$dcov <- rbind(out$dcov,bart.list[[chn]]$dcov)
##         out$dcor <- rbind(out$dcor,bart.list[[chn]]$dcor)
##         out$dsY <- rbind(out$dsY,bart.list[[chn]]$dsY)
##         if('dbeta' %in% names(out)) out$dbeta <- c(out$dbeta,bart.list[[chn]]$dbeta)
##         if('dbeta0' %in% names(out)) out$debta0 <- rbind(out$dbeta0,bart.list[[chn]]$dbeta0)
##         if('dbeta1' %in% names(out)) out$dbeta1 <- rbind(out$dbeta1,bart.list[[chn]]$dbeta1)
##         if('dh' %in% names(out)) out$dh <- rbind(out$dh,bart.list[[chn]]$dh)
##         if('dh0' %in% names(out)) out$dh2 <- rbind(out$dh0,bart.list[[chn]]$dh0)
##         if('dh1' %in% names(out)) out$dh1 <- rbind(out$dh1,bart.list[[chn]]$dh1)
##         out$df <- rbind(out$df,bart.list[[chn]]$df)        
##         out$yhat <- rbind(out$yhat,bart.list[[chn]]$yhat)
##         out$dLL <- rbind(out$dLL,bart.list[[chn]]$dLL)
##         if('dh0.test' %in% names(out)) out$dh0.test <- rbind(out$dh0.test,bart.list[[chn]]$dh0.test)
##         if('dh1.test' %in% names(out)) out$dh1.test <- rbind(out$dh1.test,bart.list[[chn]]$dh1.test)
##         out$df.test <- rbind(out$df.test,bart.list[[chn]]$df.test)
##         if('yhat.test' %in% names(out)) out$yhat.test <- rbind(out$yhat.test,bart.list[[chn]]$yhat.test)
##         out$dvarcount <- rbind(out$dvarcount,bart.list[[chn]]$dvarcount)
##         out$dvarprob <- rbind(out$dvarprob,bart.list[[chn]]$dvarprob)
##         out$chain <- c(out$chain,rep(chn,nrow(bart.list[[chn]]$dmu1)))
##     }
##     #out$yhat.mean <- colMeans(out$yhat)
##     #out$yhat.test.mean <- colMeans(out$yhat.test)    
##     return(out)
## }
        
        
