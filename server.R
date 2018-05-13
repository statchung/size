library(shiny)
library(shinythemes)
library(fpow) 

function(input, output,session) {
    observeEvent(input$do, {
    #fsize.m
      main_list<-list()
      for(i in 1:input$nf)
      {
      main_list[[i]]<-toupper(letters)[i]
      }
      
     
      
      if(max(input$checkGroup)==1)
      {
        full_list<-main_list
        updateSelectInput(session,"plot_order",choices=full_list)
      }
      else if(max(input$checkGroup)==2)
      {
        k<-1
        two_list<-list()
        
        for(i in 1:(length(main_list)-1)) {
          for(j in 1:(length(main_list)-i)){
            two_list[[k]]<-paste(main_list[[i]],"*",main_list[[i+j]])
            k<-k+1
          }
        }
        
        full_list<-c(main_list,two_list)
      
       updateSelectInput(session,"plot_order",choices= full_list)
      }
      
    fsize <-  function(alpha, beta, nu1, nu2, c,delta_type){
      fc <- qf(1-alpha, nu1, nu2)
      fl <- ncparamF(alpha, beta, nu1, nu2)/2
      if (delta_type==1) (Delta<-sqrt(2*fl/(c*nu1))) #ncp<-Delta^2*c*nu1
      else if (delta_type==2) (Delta<-sqrt(4*fl/c)) #ncp<-Delta^2*c/2
      
      return(Delta)
    }
    
    sampleSize.Factorial <- function(factor, factor.lev,delta_type, order=c(1,2),  Deltao=c(1,1), alpha=0.05, beta=0.2){
      
      main_n<-0
      two_n<-0
      nn<-0
      Delta <- NULL
      
      for (n in 2:100){
        
        v1=n-1
        if (order==1){
          v <- factor.lev-1
          c <- prod(factor.lev)*n/factor.lev
          
          v.denom <- prod(factor.lev)*(n-1)
          
          
          for (i in 1: length(v)){
            Delta[i] <- fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
          }	
          
          if (max(Delta)<=Deltao[1] ) (nn<-n)
          
        }
        else if (order==2) {
          v <- (factor.lev-1)%*%t(factor.lev-1)
          v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
          c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
          v.denom <- prod(factor.lev)*(n-1)
          
          if(main_n==0)
          {
            for (i in 1: factor){
              Delta[i] <- fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
            }	
            if (max(Delta[1:factor])<=Deltao[1] ) (main_n<-n)
          }
          
          if(two_n==0)
          {
            for (i in (factor+1): length(v)){
              Delta[i] <- fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
            }	
            if (max(Delta[(factor+1):length(v)])<=Deltao[2] ) (two_n<-n)
          }
          if(main_n >0 & two_n>0) (nn<-max(main_n,two_n))
          
        }
        
        if(nn>0) 
          {
          for(i in 1:length(v)){
            Delta[i]<-fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
          }
          break
        }
      } 
     
      return(list(n=nn, Delta=Delta))
      
    }
     
      output$Size1<-renderText({sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  delta_type=input$delta_type,order=max(input$checkGroup), Deltao=ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12)), beta=1-input$b, alpha=input$a)[[1]]}) 
    output$Size2<-renderText({sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  delta_type=input$delta_type,order=max(input$checkGroup), Deltao=ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12)), beta=1-input$b, alpha=input$a)[[2]]})
    #output$Size1<-renderText(input$delta_type)
    #output$Size2<-renderPrint(ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12)))
  #output$Size1<-FF$n
  #output$Size2<-renderText(FF$Delta)
   
   #####Graph
   output$Delta_graph<-renderPlot({
     FF2<-sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  order=max(input$checkGroup), delta_type=input$delta_type,Deltao=ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12)), beta=1-input$b, alpha=input$a)
     (n.choose <- FF2$n); 
     (Delta.choose <- data.frame(t(FF2$Delta)))
     gsize<-length(Delta.choose)
     power <- round(seq(0,1,length.out=101),3)
     
     start <- ifelse((n.choose-1)<=3, 2, n.choose-1) 
     Delta <- array(0,c(100,n.choose-start+1, ncol(Delta.choose)))
     delta.pwr <- matrix(0,n.choose-start+1, ncol(Delta.choose))
     factor.lev<-as.numeric(unlist(strsplit(input$fl,",")))
     temp_v<-list()
     temp_c<-list()
     temp_denom<-list()
     temp_n<-list()
     k<-1
     for (n in start:(n.choose)) {
       
       v1=n-1
       
       if (max(input$checkGroup)==1){
         v <- factor.lev-1
         c <- prod(factor.lev)*n/factor.lev
       } else if (max(input$checkGroup)==2) {
         v <- (factor.lev-1)%*%t(factor.lev-1)
         v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
         c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
         
       }
       
       
       v.denom <- prod(factor.lev)*(n-1)
       
       for (j in 1: length(v) ){
         delta.pwr[(n-n.choose)+nrow(delta.pwr),j]=fsize(input$a, 1-input$b, v[j], v.denom, c[j],input$delta_type);
         
         for (ind in 1: 100){
           
           if(input$a+1-power[ind]<0.9999)
             (Delta[ind,(n-n.choose)+nrow(delta.pwr),j] <- fsize(input$a, 1-power[ind], v[j], v.denom, c[j],input$delta_type))
           else (Delta[ind,(n-n.choose)+nrow(delta.pwr),j]<-NA)
            
         }	}
       temp_n[[k]]<-n
       temp_v[[k]]<-v
       temp_c[[k]]<-c
       temp_denom[[k]]<-v.denom
       k<-k+1
     }
     
     
     
     x<-seq(1:length(full_list))
     i<-x[full_list==input$plot_order]
     if(max(input$checkGroup)==1 || i<=input$nf)
     (
     de<- ifelse(input$delta_type==1,input$de1,input$de11 ) 
     )
     else if (max(input$checkGroup)==2 & i>input$nf)
       (de<-ifelse(input$delta_type==1,input$de2,input$de12 ) )
      plot(Delta[,nrow(delta.pwr),i], power[1:100], type="l", ylab="Power", xlab=paste0(input$plot_order,"_Delta"), main=paste0(input$plot_order,"_Delta vs Power(=", input$b,")"),col=1, lwd=2)
          abline(h=input$b, v=delta.pwr[nrow(delta.pwr),i], col="gray", lty=3)
         abline(h=power[which.min(abs(Delta[,nrow(delta.pwr),i]-de))], v=de, col="gray", lty=4, lwd=1.5)
         
    if(n.choose>2)
    {
         points(Delta[,nrow(delta.pwr)-1,i], power[1:100], type="l", lty=2, lwd=2)
          abline(h=input$b, v=delta.pwr[nrow(delta.pwr)-1,i], col="gray", lty=3)
          abline(h=power[which.min(abs(Delta[,nrow(delta.pwr)-1,i]-de))], v=de, col="gray", lty=4, lwd=1.5)
          
    }
         
         if(n.choose>2)
         {    
         sliderValues <- reactive({
           
           data.frame(
             R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4)),
             SD = c(round(fsize(input$a,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
                    round(fsize(input$a,0.1,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
               "1.0","1.5",
               round(fsize(input$a,0.2,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type),3),
               round(fsize(input$a,0.1,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type),3),
               "1.0","1.5"
               ),
             Power = c("0.8","0.9",
                       round((1-pf(qf((1-input$a),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type==1,(1*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                        (1*temp_c[[1]][i]/2)))),3),#(Delta^2)*(c*nu1),
                             
                       round((1-pf(qf((1-input$a),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                            (1.5^2*temp_c[[1]][i]/2)))),3),
                       "0.8","0.9",
                       round((1-pf(qf((1-input$a),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type==1,(1*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                            (1*temp_c[[2]][i]/2)))),3),#(Delta^2)*(c*nu1)
                       round((1-pf(qf((1-input$a),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type==1,(1.5^2*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                                (1.5^2*temp_c[[2]][i]/2)))),3)),
             stringsAsFactors = FALSE)
           
         })
         }
         else if (n.choose==2)
         {    
           sliderValues <- reactive({
             
             data.frame(
               R=c(rep(temp_n[[1]],4)),
               SD = c(round(fsize(input$a,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
                      round(fsize(input$a,0.1,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
                      "1.0","1.5"
               ),
               Power = c("0.8","0.9",
                         round((1-pf(qf((1-input$a),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type==1,(1*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                              (1*temp_c[[1]][i]/2)))),3),#(Delta^2)*(c*nu1)
                         round((1-pf(qf((1-input$a),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                              (1.5^2*temp_c[[1]][i]/2)))),3)
                         ),
               stringsAsFactors = FALSE)
             
           })
         }
         # Show the values in an HTML table ----
         output$values <- renderTable({
           sliderValues()
         })
         
         
   })
   
   
   
   
   output$Size_graph<-renderPlot({
     FF2<-sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  order=max(input$checkGroup), delta_type=input$delta_type,Deltao=ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12)), beta=1-input$b, alpha=input$a)
 
     (n.choose <- FF2$n); 
     (Delta.choose <- data.frame(t(FF2$Delta)))
     gsize<-length(Delta.choose)
     power <- round(seq(0,1,length.out=101),3)
 
     factor.lev<-as.numeric(unlist(strsplit(input$fl,",")))
     
     Delta <- matrix(0,100,ncol(Delta.choose)) 
     
     for (n in 2:100) {
       
       v1=n-1
      
       if (max(input$checkGroup)==1){
         v <- factor.lev-1
         c <- prod(factor.lev)*n/factor.lev
       } else if (max(input$checkGroup)==2) {
         v <- (factor.lev-1)%*%t(factor.lev-1)
         v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
         c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
         
       }
       
       
       v.denom <- prod(factor.lev)*(n-1)
       
       for (j in 1: length(v) ){
         Delta[n,j] <- fsize(input$a, 1-input$b, v[j], v.denom, c[j],input$delta_type)
         	}
     }
     
 
         plot(2:100, Delta[2:100,1], type="l", xlim=c(0,min(50,n.choose+5)), ylim=c(0,max(ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12)))*1.5), ylab="Delta", xlab="the number of replications", 
            main="Delta vs the number of replications(r)",col=1, lwd=2)
       for (i in 2:ncol(Delta)) 
         lines(2:100, Delta[2:100,i], type="l", lty=i, lwd=2,col=i)
          
       abline(h=max(ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12))), v=FF2$n,col="gray", lty=3)
       legend("top", legend=paste0("power=", input$b), adj=0, bty="n")
       
       legend("topright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
       
       
   })
   
   
   output$power_graph<-renderPlot({
     FF2<-sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  order=max(input$checkGroup), delta_type=input$delta_type,Deltao=ifelse(rep(input$delta_type==1,2),c(input$de1,input$de2),c(input$de11,input$de12)), beta=1-input$b, alpha=input$a)
     (n.choose <- FF2$n); 
     (Delta.choose <- data.frame(t(FF2$Delta)))
     gsize<-length(Delta.choose)
     power <- round(seq(0,1,length.out=101),3)
     Deltao <- c(1,1.5,2)
     Delta <- pwr <- array(0,c(100, ncol(Delta.choose),3))
     delta.pwr <- array(0,c(100,ncol(Delta.choose),3))
     factor.lev<-as.numeric(unlist(strsplit(input$fl,",")))
     
     
     for (deltao in 1: 3){
       for (n in 2:100) {
       
       v1=n-1
       
       if (max(input$checkGroup)==1){
         v <- factor.lev-1
         c <- prod(factor.lev)*n/factor.lev
       } else if (max(input$checkGroup)==2) {
         v <- (factor.lev-1)%*%t(factor.lev-1)
         v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
         c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
         
       }
       
       
       v.denom <- prod(factor.lev)*(n-1)
         
         for (j in 1: length(v)){ 
          pwr[n,j,deltao]<-round((1-pf(qf((1-input$a),v[j],v.denom),v[j],v.denom,ncp=ifelse(input$delta_type==1,(Deltao[deltao]^2*(c[j]*v[j])),
                                                                                            (Deltao[deltao]^2*c[j]/2)))),3)
                                       
                                        
         } 
    
       
       
       
       }
     }
      
     x<-c(1,2,3)
     j<-x[Deltao==input$plot_delta]
   
       plot(2:100, pwr[2:100,1, j], type="l", ylim=c(0,1), xlim=c(0,min(50,n.choose+5)), ylab="power", xlab="the number of replications", 
            main="power vs the number of replications(r) ",col=1, lwd=2)
       for (i in 2:ncol(Delta)) 
         lines(2:100, pwr[2:100,i,j], type="l", lty=i, lwd=2,col=i)
       
       abline(h=0.8, v=FF2$n,col="gray", lty=3)
       
       abline(h=0.9, col="gray", lty=3)
       legend("bottomright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
       
       }
     
   
   )
    }
 )
  }
