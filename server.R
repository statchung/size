library(shiny)
library(shinythemes)
library(fpow) 


  server <- function(input, output,session) {
    observeEvent(input$do, {
    #fsize.m
      main_list<-list()
      for(i in 1:input$nf)
      {
      main_list[[i]]<-toupper(letters)[i]
      }
      
     
      
      if(max(input$checkGroup)==1)
      {
        full_list1<-c(main_list,"ALL")
        full_list<-main_list
        updateSelectInput(session,"plot_order",choices=full_list1)
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
        
        full_list1<-c(main_list,two_list,"ALL")
        full_list<-c(main_list,two_list)
       updateSelectInput(session,"plot_order",choices= full_list1)
      }
      
    fsize <-  function(alpha, beta, nu1, nu2, c,delta_type){
      fc <- qf(1-alpha, nu1, nu2)
      fl <- ncparamF(alpha, beta, nu1, nu2)/2
      if (delta_type==1) (Delta<-sqrt(2*fl/(c*nu1))) #ncp<-Delta^2*c*nu1
      else if (delta_type==2) (Delta<-sqrt(4*fl/c)) #ncp<-Delta^2*c/2
      
      return(Delta)
    }
    
    list_fc<-function(list)
    {
      x<- length(list) 
      
      list1<-list[[1]]
      for( i in 1:(x-1))
        list1<-paste(list1,list[[i+1]],sep="+")
      return(list1=list1)
    }
    
    sampleSize.Factorial <- function(factor, factor.lev,delta_type, order=c(1,2),  Deltao=c(1,1,1), alpha=0.05, beta=0.2){
      
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
          
          if (max(Delta)<=Deltao[1]/Deltao[3] ) (nn<-n)
          
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
            if (max(Delta[1:factor])<=Deltao[1]/Deltao[3] ) (main_n<-n)
          }
          
          if(two_n==0)
          {
            for (i in (factor+1): length(v)){
              Delta[i] <- fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
            }	
            if (max(Delta[(factor+1):length(v)])<=Deltao[2]/Deltao[3] ) (two_n<-n)
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
   
     
    output$Size1<-renderText({sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  delta_type=input$delta_type,order=max(input$checkGroup), Deltao=ifelse(rep(input$delta_type==1,3),c(input$de1,input$de2,input$de3),c(input$de11,input$de12,input$de13)), beta=1-input$b, alpha=input$a)[[1]]}) 
    output$Size2<-renderText({sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  delta_type=input$delta_type,order=max(input$checkGroup), Deltao=ifelse(rep(input$delta_type==1,3),c(input$de1,input$de2,input$de3),c(input$de11,input$de12,input$de13)), beta=1-input$b, alpha=input$a)[[2]]})
    output$list1<-renderText({list_fc(full_list)})
    
   #####Graph
   output$Delta_graph<-renderPlot({
     FF2<-sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  order=max(input$checkGroup), delta_type=input$delta_type,Deltao=ifelse(rep(input$delta_type==1,3),c(input$de1,input$de2,input$de3),c(input$de11,input$de12,input$de13)), beta=1-input$b, alpha=input$a)
     (n.choose <- FF2$n); 
     (Delta.choose <- data.frame(t(FF2$Delta)))
     gsize<-length(Delta.choose)
     power <- round(seq(0,1,length.out=101),3)
     
     start <- ifelse((n.choose-1)<=2, 2, n.choose-2) 
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
     
     
     
     x<-seq(1:length(full_list1))
     i<-x[full_list1==input$plot_order] 
     
     if(max(input$checkGroup)==1 || i<=input$nf)
     (
     de<- ifelse(input$delta_type==1,input$de1/input$de3,input$de11/input$de13 ) 
     )
     else if (max(input$checkGroup)==2 & i>input$nf)
       (de<-ifelse(input$delta_type==1,input$de2/input$de3,input$de12/input$de13 ) )
     if(i==max(x))
     {
       plot(Delta[,nrow(delta.pwr),1], power[1:100],  type="l", ylab="Power", 
            xlab=ifelse(input$plot_order=="ALL","Delta",ifelse(input$delta_type==1,paste0("SD(",input$plot_order,")/SD(noise)"),paste0("Range(",input$plot_order,")/SD(noise)"))) ,
            main=ifelse(input$plot_order=="ALL",paste0("Delta vs Power(=",input$b,")"),ifelse(input$delta_type==1,paste0("SD(",input$plot_order,")/SD(noise) vs Power(=", input$b,")"),paste0("Range(",input$plot_order,")/SD(noise) vs Power(=", input$b,")")))
            ,col=1, lty=1,lwd=2)
       for(j in 2:(i-1))
       {
       
         points(Delta[,nrow(delta.pwr),j], power[1:100], type="l", col=j,lty=j, lwd=2)
       }
       abline(h=c(0.8,0.9),v=c(1.0,1.5),col="gray")
       
       legend("bottomleft", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
       
       sliderValues <- reactive({})
       
       output$values <- renderTable({
         sliderValues()
       })
     }
     
     
     
     else if(i<max(x))
     {     
      plot(Delta[,nrow(delta.pwr),i], power[1:100], type="l", ylab="Power", 
           xlab=ifelse(input$plot_order=="ALL","Delta",ifelse(input$delta_type==1,paste0("SD(",input$plot_order,")/SD(noise)"),paste0("Range(",input$plot_order,")/SD(noise)"))) ,
           main=ifelse(input$plot_order=="ALL",paste0("Delta vs Power(=",input$b,")"),ifelse(input$delta_type==1,paste0("SD(",input$plot_order,")/SD(noise) vs Power(=", input$b,")"),paste0("Range(",input$plot_order,")/SD(noise) vs Power(=", input$b,")")))
           ,col=1, lwd=2)
         abline(h=0.8, v= fsize(input$a,0.2,temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],temp_c[[nrow(delta.pwr)]][i],input$delta_type) , col="gray", lty=3)
          abline(h=0.9, v =fsize(input$a,0.1,temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],temp_c[[nrow(delta.pwr)]][i],input$delta_type) , col="gray", lty=3)
          abline(h= (1-pf(qf((1-input$a),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]]),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],ncp=ifelse(input$delta_type==1,(1*(temp_c[[nrow(delta.pwr)]][i]*temp_v[[nrow(delta.pwr)]][i])),
                                                                                                                        (1*temp_c[[nrow(delta.pwr)]][i]/2)))) , v =1.0, col="gray", lty=3)
          abline(h= (1-pf(qf((1-input$a),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]]),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],ncp=ifelse(input$delta_type==1,(1.5^2*(temp_c[[nrow(delta.pwr)]][i]*temp_v[[nrow(delta.pwr)]][i])),
                                                                                                                        (1.5^2*temp_c[[nrow(delta.pwr)]][i]/2)))) , v =1.5 , col="gray", lty=3)
         
         
    if(n.choose>3)
    {
         points(Delta[,nrow(delta.pwr)-1,i], power[1:100], type="l", lty=2, lwd=2)
      
      points(Delta[,nrow(delta.pwr)-2,i], power[1:100], type="l", lty=3, lwd=2)
       
          
    }
          
         else if(n.choose==3)
          {
            points(Delta[,nrow(delta.pwr)-1,i], power[1:100], type="l", lty=2, lwd=2) 
            
            
          }
     
         if(n.choose>3)
         {    
         sliderValues <- reactive({
           
           data.frame(
             R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4),rep(temp_n[[3]],4)),
             Delta = c(round(fsize(input$a,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
                    round(fsize(input$a,0.1,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
               "1.0","1.5",
               round(fsize(input$a,0.2,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type),3),
               round(fsize(input$a,0.1,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type),3),
               "1.0","1.5",
               
               round(fsize(input$a,0.2,temp_v[[3]][i],temp_denom[[3]],temp_c[[3]][i],input$delta_type),3),
               round(fsize(input$a,0.1,temp_v[[3]][i],temp_denom[[3]],temp_c[[3]][i],input$delta_type),3),
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
                                                                                                                                (1.5^2*temp_c[[2]][i]/2)))),3),
             "0.8","0.9",
             round((1-pf(qf((1-input$a),temp_v[[3]][i],temp_denom[[3]]),temp_v[[3]][i],temp_denom[[3]],ncp=ifelse(input$delta_type==1,(1*(temp_c[[3]][i]*temp_v[[3]][i])),
                                                                                                                  (1*temp_c[[3]][i]/2)))),3),#(Delta^2)*(c*nu1)
             round((1-pf(qf((1-input$a),temp_v[[3]][i],temp_denom[[3]]),temp_v[[3]][i],temp_denom[[3]],ncp=ifelse(input$delta_type==1,(1.5^2*(temp_c[[3]][i]*temp_v[[3]][i])),
                                                                                                                  (1.5^2*temp_c[[3]][i]/2)))),3)),
             stringsAsFactors = FALSE)
           
         })
         }
         else if (n.choose==3)
         {    
           sliderValues <- reactive({
             
             data.frame(
               R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4)),
               Delta = c(round(fsize(input$a,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
                      round(fsize(input$a,0.1,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
                      "1.0","1.5",
                      round(fsize(input$a,0.2,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type),3),
                      round(fsize(input$a,0.1,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type),3),
                      "1.0","1.5"
               ),
               Power = c("0.8","0.9",
                         round((1-pf(qf((1-input$a),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type==1,(1*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                              (1*temp_c[[1]][i]/2)))),3),#(Delta^2)*(c*nu1)
                         round((1-pf(qf((1-input$a),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                              (1.5^2*temp_c[[1]][i]/2)))),3),
                         "0.8","0.9",
                         round((1-pf(qf((1-input$a),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type==1,(1*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                              (1*temp_c[[2]][i]/2)))),3),#(Delta^2)*(c*nu1)
                         round((1-pf(qf((1-input$a),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type==1,(1.5^2*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                              (1.5^2*temp_c[[2]][i]/2)))),3)
                         
                         ),
               stringsAsFactors = FALSE)
             
           })
         }
          else if (n.choose==2)
          {    
            sliderValues <- reactive({
              
              data.frame(
                R=c(rep(temp_n[[1]],4)),
                Delta = c(round(fsize(input$a,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type),3),
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
     }
         
   })
   
   
   
   
   output$Size_graph<-renderPlot({
     FF2<-sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  order=max(input$checkGroup), delta_type=input$delta_type,Deltao=ifelse(rep(input$delta_type==1,3),c(input$de1,input$de2,input$de3),c(input$de11,input$de12,input$de13)), beta=1-input$b, alpha=input$a)
 
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
     
 
         plot(2:100, Delta[2:100,1], type="l", xlim=c(0,min(100,n.choose+5)), ylim=c(0,max(ifelse(rep(input$delta_type==1,2),c(input$de1/input$de3,input$de2/input$de3),c(input$de11/input$de13,input$de12/input$de13)))*1.5), ylab="Delta", xlab="The number of replications", 
            main="The number of replications(r) vs Delta",col=1, lwd=2)
       for (i in 2:ncol(Delta)) 
         lines(2:100, Delta[2:100,i], type="l", lty=i, lwd=2,col=i)
          
       abline(h=max(ifelse(rep(input$delta_type==1,2),c(input$de1/input$de3,input$de2/input$de3),c(input$de11/input$de13,input$de12/input$de13))), v=FF2$n,col="gray", lty=3)
       legend("top", legend=paste0("power=", input$b), adj=0, bty="n")
       
       legend("topright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
       
       
   })
   
   
   output$power_graph<-renderPlot({
     FF2<-sampleSize.Factorial(input$nf, as.numeric(unlist(strsplit(input$fl,","))),  order=max(input$checkGroup), delta_type=input$delta_type,Deltao=ifelse(rep(input$delta_type==1,3),c(input$de1,input$de2,input$de3),c(input$de11,input$de12,input$de13)), beta=1-input$b, alpha=input$a)
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
   
       plot(2:100, pwr[2:100,1, j], type="l", ylim=c(0,1), xlim=c(0,min(100,n.choose+5)), ylab="Power", xlab="The number of replications", 
            main="The number of replications(r) vs Power",col=1, lwd=2)
       for (i in 2:ncol(Delta)) 
         lines(2:100, pwr[2:100,i,j], type="l", lty=i, lwd=2,col=i)
       
       abline(h=0.8, v=FF2$n,col="gray", lty=3)
       
       abline(h=0.9, col="gray", lty=3)
       legend("bottomright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
       
       }
     
   
   )
    }
 ) 
    
    observeEvent(input$do2, {
      #fsize.m
     
      main_list<-list()
      for(i in 1:input$nf2)
      {
        main_list[[i]]<-toupper(letters)[i]
      }
      
      
      fsize <-  function(alpha, beta, nu1, nu2, c,delta_type){
        fc <- qf(1-alpha, nu1, nu2)
        fl <- ncparamF(alpha, beta, nu1, nu2)/2
        if (delta_type==1) (Delta<-sqrt(2*fl/(c*nu1))) #ncp<-Delta^2*c*nu1
        else if (delta_type==2) (Delta<-sqrt(4*fl/c)) #ncp<-Delta^2*c/2
        
        return(Delta)
      }
      
      list_fc<-function(list)
      {
        x<- length(list) 
        
        list1<-list[[1]]
        for( i in 1:(x-1))
          list1<-paste(list1,list[[i+1]],sep="+")
        return(list1=list1)
      }
      
      sampleSize.2levFrFactorial <- function(factor, fraction, resolution, delta_type, Deltao=c(1,1), alpha=0.05, beta=0.2){
        
        for (n in 2:100){
          
          v1=n-1
          
          v2=1
          v.denom <- ifelse(resolution==3, 2^(factor-fraction)*n-1-1, 2^(factor-fraction)*n-1-factor)
          
          Ac=2^(factor-fraction-1)*n
          
          ADelta=fsize(alpha, beta, v2, v.denom, Ac,delta_type);
          
          if (ADelta<=Deltao[1]/Deltao[2] ) break
          
        }
        
        return(list(n=n, Delta=ADelta))
        
      }
      
      
      
      output$Size12<-renderText({sampleSize.2levFrFactorial(input$nf2, input$fr2, input$check_res,delta_type=input$delta_type2,Deltao=ifelse(rep(input$delta_type2==1,2),c(input$de1_2,input$de3_2),c(input$de11_2,input$de13_2)), beta=1-input$b2, alpha=input$a2)[[1]]})
      output$Size22<-renderText({sampleSize.2levFrFactorial(input$nf2, input$fr2, input$check_res,delta_type=input$delta_type2,Deltao=ifelse(rep(input$delta_type2==1,2),c(input$de1_2,input$de3_2),c(input$de11_2,input$de13_2)), beta=1-input$b2, alpha=input$a2)[[2]]})
      output$list1_2<-renderText({list_fc(main_list)})
      
      
      output$Delta_graph2<-renderPlot({
        FF2<-sampleSize.2levFrFactorial(input$nf2, input$fr2, input$check_res,delta_type=input$delta_type2,Deltao=ifelse(rep(input$delta_type2==1,2),c(input$de1_2,input$de3_2),c(input$de11_2,input$de13_2)), beta=1-input$b2, alpha=input$a2)
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        
        start <- ifelse((n.choose-1)<=2, 2, n.choose-2) 
        Delta <- array(0,c(100,n.choose-start+1, ncol(Delta.choose)))
        delta.pwr <- matrix(0,n.choose-start+1, ncol(Delta.choose))
        factor.lev<-2
        temp_v<-list()
        temp_c<-list()
        temp_denom<-list()
        temp_n<-list()
        k<-1
        
        for (n in start:(n.choose)) {
          
          v1=n-1
          
          v=1
          v.denom <- ifelse(input$check_res==3, 2^(input$nf2-input$fr2)*n-1-1, 2^(input$nf2-input$fr2)*n-1-input$nf2)
          
          c=2^(input$nf2-input$fr2-1)*n
          
          delta.pwr[(n-n.choose)+nrow(delta.pwr),]=fsize(input$a2, 1-input$b2, v, v.denom, c,input$delta_type2);
          
          for (ind in 1: 100){
            if(input$a2+1-power[ind]<0.9999)
              (Delta[ind,(n-n.choose)+nrow(delta.pwr),] <- fsize(input$a2, 1-power[ind], v, v.denom, c,input$delta_type2))
            else (Delta[ind,(n-n.choose)+nrow(delta.pwr),]<-NA)
          }	
          temp_n[[k]]<-n
          temp_v[[k]]<-v
          temp_c[[k]]<-c
          temp_denom[[k]]<-v.denom
          k<-k+1
        }
        
        
        
        
        de<-ifelse(input$delta_type2==1,input$de1_2/input$de3_2,input$de11_2/input$de13_2)
        
        
        plot(Delta[,nrow(delta.pwr),1], power[1:100], type="l", ylab="Power", xlab=ifelse(input$delta_type2==1,"SD(Main)/SD(noise)", "Range(Main)/SD(noise)" ) ,
             main=ifelse(input$delta_type2==1,paste0("SD(Main)/SD(noise) vs Power(=", input$b2,")"),paste0("Range(Main)/SD(noise) vs Power(=", input$b2,")"))
             ,col=1, lty=1,lwd=2)
        
        abline(h=0.8, v= fsize(input$a2,0.2,temp_v[[nrow(delta.pwr)]],temp_denom[[nrow(delta.pwr)]],temp_c[[nrow(delta.pwr)]],input$delta_type2) , col="gray", lty=3)
        abline(h=0.9, v =fsize(input$a2,0.1,temp_v[[nrow(delta.pwr)]],temp_denom[[nrow(delta.pwr)]],temp_c[[nrow(delta.pwr)]],input$delta_type2) , col="gray", lty=3)
        abline(h= (1-pf(qf((1-input$a2),temp_v[[nrow(delta.pwr)]],temp_denom[[nrow(delta.pwr)]]),temp_v[[nrow(delta.pwr)]],temp_denom[[nrow(delta.pwr)]],ncp=ifelse(input$delta_type2==1,(1*(temp_c[[nrow(delta.pwr)]]*temp_v[[nrow(delta.pwr)]])),
                                                                                                                                                                   (1*temp_c[[nrow(delta.pwr)]]/2)))) , v =1.0, col="gray", lty=3)
        abline(h= (1-pf(qf((1-input$a2),temp_v[[nrow(delta.pwr)]],temp_denom[[nrow(delta.pwr)]]),temp_v[[nrow(delta.pwr)]],temp_denom[[nrow(delta.pwr)]],ncp=ifelse(input$delta_type2==1,(1.5^2*(temp_c[[nrow(delta.pwr)]]*temp_v[[nrow(delta.pwr)]])),
                                                                                                                                                                   (1.5^2*temp_c[[nrow(delta.pwr)]]/2)))) , v =1.5 , col="gray", lty=3)
        
        if(n.choose>3)
        {
          points(Delta[,nrow(delta.pwr)-1,1], power[1:100], type="l", lty=2, lwd=2)
          
          points(Delta[,nrow(delta.pwr)-2,1], power[1:100], type="l", lty=3, lwd=2)
          
          
        }
        
        else if(n.choose==3)
        {
          points(Delta[,nrow(delta.pwr)-1,1], power[1:100], type="l", lty=2, lwd=2) 
          
          
        }
        
        if(n.choose>3)
        {    
          sliderValues2 <- reactive({
            
            data.frame(
              R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4),rep(temp_n[[3]],4)),
              Delta = c(round(fsize(input$a2,0.2,temp_v[[1]],temp_denom[[1]],temp_c[[1]],input$delta_type2),3),
                        round(fsize(input$a2,0.1,temp_v[[1]],temp_denom[[1]],temp_c[[1]],input$delta_type2),3),
                        "1.0","1.5",
                        round(fsize(input$a2,0.2,temp_v[[2]],temp_denom[[2]],temp_c[[2]],input$delta_type2),3),
                        round(fsize(input$a2,0.1,temp_v[[2]],temp_denom[[2]],temp_c[[2]],input$delta_type2),3),
                        "1.0","1.5",
                        
                        round(fsize(input$a2,0.2,temp_v[[3]],temp_denom[[3]],temp_c[[3]],input$delta_type2),3),
                        round(fsize(input$a2,0.1,temp_v[[3]],temp_denom[[3]],temp_c[[3]],input$delta_type2),3),
                        "1.0","1.5"
              ),
              Power = c("0.8","0.9",
                        round((1-pf(qf((1-input$a2),temp_v[[1]],temp_denom[[1]]),temp_v[[1]],temp_denom[[1]],ncp=ifelse(input$delta_type2==1,(1*(temp_c[[1]]*temp_v[[1]])),
                                                                                                                       (1*temp_c[[1]]/2)))),3),#(Delta^2)*(c*nu1),
                        
                        round((1-pf(qf((1-input$a2),temp_v[[1]],temp_denom[[1]]),temp_v[[1]],temp_denom[[1]],ncp=ifelse(input$delta_type2==1,(1.5^2*(temp_c[[1]]*temp_v[[1]])),
                                                                                                                       (1.5^2*temp_c[[1]]/2)))),3),
                        "0.8","0.9",
                        round((1-pf(qf((1-input$a2),temp_v[[2]],temp_denom[[2]]),temp_v[[2]],temp_denom[[2]],ncp=ifelse(input$delta_type2==1,(1*(temp_c[[2]]*temp_v[[2]])),
                                                                                                                       (1*temp_c[[2]]/2)))),3),#(Delta^2)*(c*nu1)
                        round((1-pf(qf((1-input$a2),temp_v[[2]],temp_denom[[2]]),temp_v[[2]],temp_denom[[2]],ncp=ifelse(input$delta_type2==1,(1.5^2*(temp_c[[2]]*temp_v[[2]])),
                                                                                                                       (1.5^2*temp_c[[2]]/2)))),3),
                        "0.8","0.9",
                        round((1-pf(qf((1-input$a2),temp_v[[3]],temp_denom[[3]]),temp_v[[3]],temp_denom[[3]],ncp=ifelse(input$delta_type2==1,(1*(temp_c[[3]]*temp_v[[3]])),
                                                                                                                       (1*temp_c[[3]]/2)))),3),#(Delta^2)*(c*nu1)
                        round((1-pf(qf((1-input$a2),temp_v[[3]],temp_denom[[3]]),temp_v[[3]],temp_denom[[3]],ncp=ifelse(input$delta_type2==1,(1.5^2*(temp_c[[3]]*temp_v[[3]])),
                                                                                                                       (1.5^2*temp_c[[3]]/2)))),3)),
              stringsAsFactors = FALSE)
            
          })
        }
        else if (n.choose==3)
        {    
          sliderValues2 <- reactive({
            
            data.frame(
              R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4)),
              Delta = c(round(fsize(input$a2,0.2,temp_v[[1]],temp_denom[[1]],temp_c[[1]],input$delta_type2),3),
                        round(fsize(input$a2,0.1,temp_v[[1]],temp_denom[[1]],temp_c[[1]],input$delta_type2),3),
                        "1.0","1.5",
                        round(fsize(input$a2,0.2,temp_v[[2]],temp_denom[[2]],temp_c[[2]],input$delta_type2),3),
                        round(fsize(input$a2,0.1,temp_v[[2]],temp_denom[[2]],temp_c[[2]],input$delta_type2),3),
                        "1.0","1.5"
              ),
              Power = c("0.8","0.9",
                        round((1-pf(qf((1-input$a2),temp_v[[1]],temp_denom[[1]]),temp_v[[1]],temp_denom[[1]],ncp=ifelse(input$delta_type2==1,(1*(temp_c[[1]]*temp_v[[1]])),
                                                                                                                       (1*temp_c[[1]]/2)))),3),#(Delta^2)*(c*nu1)
                        round((1-pf(qf((1-input$a2),temp_v[[1]],temp_denom[[1]]),temp_v[[1]],temp_denom[[1]],ncp=ifelse(input$delta_type2==1,(1.5^2*(temp_c[[1]]*temp_v[[1]])),
                                                                                                                       (1.5^2*temp_c[[1]]/2)))),3),
                        "0.8","0.9",
                        round((1-pf(qf((1-input$a2),temp_v[[2]],temp_denom[[2]]),temp_v[[2]],temp_denom[[2]],ncp=ifelse(input$delta_type2==1,(1*(temp_c[[2]]*temp_v[[2]])),
                                                                                                                       (1*temp_c[[2]]/2)))),3),#(Delta^2)*(c*nu1)
                        round((1-pf(qf((1-input$a2),temp_v[[2]],temp_denom[[2]]),temp_v[[2]],temp_denom[[2]],ncp=ifelse(input$delta_type2==1,(1.5^2*(temp_c[[2]]*temp_v[[2]])),
                                                                                                                       (1.5^2*temp_c[[2]]/2)))),3)
                        
              ),
              stringsAsFactors = FALSE)
            
          }) 
          
         
        }
        else if (n.choose==2)
        {    
          sliderValues2 <- reactive({
            
            data.frame(
              R=c(rep(temp_n[[1]],4)),
              Delta = c(round(fsize(input$a2,0.2,temp_v[[1]],temp_denom[[1]],temp_c[[1]],input$delta_type2),3),
                        round(fsize(input$a2,0.1,temp_v[[1]],temp_denom[[1]],temp_c[[1]],input$delta_type2),3),
                        "1.0","1.5" 
              ),
              Power = c("0.8","0.9",
                        round((1-pf(qf((1-input$a2),temp_v[[1]],temp_denom[[1]]),temp_v[[1]],temp_denom[[1]],ncp=ifelse(input$delta_type2==1,(1*(temp_c[[1]]*temp_v[[1]])),
                                                                                                                        (1*temp_c[[1]]/2)))),3),#(Delta^2)*(c*nu1)
                        round((1-pf(qf((1-input$a2),temp_v[[1]],temp_denom[[1]]),temp_v[[1]],temp_denom[[1]],ncp=ifelse(input$delta_type2==1,(1.5^2*(temp_c[[1]]*temp_v[[1]])),
                                                                                                                        (1.5^2*temp_c[[1]]/2)))),3) 
                        
              ),
              stringsAsFactors = FALSE)
            
          }) 
          
          
        }
        output$values2 <- renderTable({
          sliderValues2()
        })
        
        
      })
      
      output$Size_graph2<-renderPlot({
        FF2<-sampleSize.2levFrFactorial(input$nf2, input$fr2, input$check_res,delta_type=input$delta_type2,Deltao=ifelse(rep(input$delta_type2==1,2),c(input$de1_2,input$de3_2),c(input$de11_2,input$de13_2)), beta=1-input$b2, alpha=input$a2)
        
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        
        factor.lev<-2
        
        Delta <- matrix(0,100,ncol(Delta.choose)) 
        
        for (n in 2:100) {
          
          
          v1=n-1
          
          v=1
          v.denom <- ifelse(input$fr2==3, 2^(input$nf2-input$fr2)*n-1-1, 2^(input$nf2-input$fr2)*n-1-input$nf2)
          
          c=2^(input$nf2-input$fr2-1)*n
          
          Delta[n,] <- fsize(input$a2, 1-input$b2, v, v.denom, c,input$delta_type2)
          
          
        }
        
        
        plot(2:100, Delta[2:100,1], type="l", xlim=c(0,min(100,n.choose+5)), ylim=c(0,max(ifelse(input$delta_type2==1,input$de1_2/input$de3_2,input$de11_2/input$de13_2))*1.5), ylab="Delta", xlab="The number of replications", 
             main="The number of replications(r) vs Delta",col=1, lwd=2)
        
        abline(h=max(ifelse(input$delta_type2==1,input$de1_2/input$de3_2,input$de11_2/input$de13_2)), v=FF2$n,col="gray", lty=3)
        legend("top", legend=paste0("power=", input$b2), adj=0, bty="n")
        
        
      })
      
      output$power_graph2<-renderPlot({
        FF2<-sampleSize.2levFrFactorial(input$nf2, input$fr2, input$check_res,delta_type=input$delta_type2,Deltao=ifelse(rep(input$delta_type2==1,2),c(input$de1_2,input$de3_2),c(input$de11_2,input$de13_2)), beta=1-input$b2, alpha=input$a2)
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        Deltao <- c(1,1.5,2)
        Delta <- pwr <- array(0,c(100, ncol(Delta.choose),3))
        delta.pwr <- array(0,c(100,ncol(Delta.choose),3))
        factor.lev<-2
        
        
        for (deltao in 1: 3){
          for (n in 2:100) {
            
            v1=n-1
            
            v=1
            v.denom <-  ifelse(input$fr2==3, 2^(input$nf2-input$fr2)*n-1-1, 2^(input$nf2-input$fr2)*n-1-input$nf2)
            
            c=2^(input$nf2-input$fr2-1)*n
            #Delta[n,] <- fsize(alpha, beta, v, v.denom, c)
            
            pwr[n,1,deltao]<-round((1-pf(qf((1-input$a2),v ,v.denom),v ,v.denom,ncp=ifelse(input$delta_type2==1,(Deltao[deltao]^2*(c *v)),
                                                                                           (Deltao[deltao]^2*c/2)))),3)
            
            
          }
        }
        
        x<-c(1,2,3)
        j<-x[Deltao==input$plot_delta2]
        
        plot(2:100, pwr[2:100,1, j], type="l", ylim=c(0,1), xlim=c(0,min(100,n.choose+5)), ylab="Power", xlab="The number of replications", 
             main="The number of replications(r) vs Power",col=1, lwd=2)
        
        abline(h=0.8,col="gray", lty=3)
        
        abline(h=0.9, col="gray", lty=3)
        
      }
      
      
      )
    })
    
    observeEvent(input$do3, {
      #fsize.m
      main_list<-list()
      for(i in 1:input$nf3)
      {
        main_list[[i]]<-toupper(letters)[i]
      }
      
      
      
      if(max(input$checkGroup3)==1)
      {
        full_list1<-c(main_list,"ALL")
        full_list<-main_list
        updateSelectInput(session,"plot_order3",choices=full_list1)
      }
      else if(max(input$checkGroup3)==2)
      {
        k<-1
        two_list<-list()
        
        for(i in 1:(length(main_list)-1)) {
          for(j in 1:(length(main_list)-i)){
            two_list[[k]]<-paste(main_list[[i]],"*",main_list[[i+j]])
            k<-k+1
          }
        }
        
        full_list1<-c(main_list,two_list,"ALL")
        full_list<-c(main_list,two_list)
        updateSelectInput(session,"plot_order3",choices= full_list1)
      }
      
      fsize <-  function(alpha, beta, nu1, nu2, c,delta_type){
        fc <- qf(1-alpha, nu1, nu2)
        fl <- ncparamF(alpha, beta, nu1, nu2)/2
        if (delta_type==1) (Delta<-sqrt(2*fl/(c*nu1))) #ncp<-Delta^2*c*nu1
        else if (delta_type==2) (Delta<-sqrt(4*fl/c)) #ncp<-Delta^2*c/2
        
        return(Delta)
      }
      
      list_fc<-function(list)
      {
        x<- length(list) 
        
        list1<-list[[1]]
        for( i in 1:(x-1))
          list1<-paste(list1,list[[i+1]],sep="+")
        return(list1=list1)
      }
      
      sampleSize.RCBD <- function(factor, factor.lev,delta_type, order=c(1,2),  Deltao=c(1,1,1), alpha=0.05, beta=0.2){
        
        main_n<-0
        two_n<-0
        nn<-0
        Delta <- NULL
        
        for (n in 2:100){
          
          v1=n-1
          if (order==1){
            v <- factor.lev-1
            c <- prod(factor.lev)*n/factor.lev
            
            v.denom <- prod(factor.lev)*n-1-sum(v)-v1
            
            
            for (i in 1: length(v)){
              Delta[i] <- fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
            }	
            
            if (max(Delta)<=Deltao[1]/Deltao[3] ) 
            (nn<-n
            )
            
            
          }
          else if (order==2) {
            v <- (factor.lev-1)%*%t(factor.lev-1)
            v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
            c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
            v.denom <- prod(factor.lev)*n-1-sum(v)-v1
            
            if(main_n==0)
            {
              for (i in 1: factor){
                Delta[i] <- fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
              }	
              if (max(Delta[1:factor])<=Deltao[1]/Deltao[3] ) (main_n<-n)
            }
            
            if(two_n==0)
            {
              for (i in (factor+1): length(v)){
                Delta[i] <- fsize(alpha, beta, v[i], v.denom, c[i],delta_type)
              }	
              if (max(Delta[(factor+1):length(v)])<=Deltao[2]/Deltao[3] ) (two_n<-n)
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
      
      output$Size1_3<-renderText({sampleSize.RCBD(input$nf3, as.numeric(unlist(strsplit(input$fl3,","))),  delta_type=input$delta_type3,order=max(input$checkGroup3), Deltao=ifelse(rep(input$delta_type3==1,3),c(input$de1_3,input$de2_3,input$de3_3),c(input$de11_3,input$de12_3,input$de13_3)), beta=1-input$b3, alpha=input$a3)[[1]]}) 
      output$Size2_3<-renderText({sampleSize.RCBD(input$nf3, as.numeric(unlist(strsplit(input$fl3,","))),  delta_type=input$delta_type3,order=max(input$checkGroup3), Deltao=ifelse(rep(input$delta_type3==1,3),c(input$de1_3,input$de2_3,input$de3_3),c(input$de11_3,input$de12_3,input$de13_3)), beta=1-input$b3, alpha=input$a3)[[2]]})
      output$list1_3<-renderText({list_fc(full_list)})
      
      #####Graph
      output$Delta_graph3<-renderPlot({
        FF2<-sampleSize.RCBD(input$nf3, as.numeric(unlist(strsplit(input$fl3,","))),  order=max(input$checkGroup3), delta_type=input$delta_type3,Deltao=ifelse(rep(input$delta_type3==1,3),c(input$de1_3,input$de2_3,input$de3_3),c(input$de11_3,input$de12_3,input$de13_3)), beta=1-input$b3, alpha=input$a3)
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        
        start <- ifelse((n.choose-1)<=2, 2, n.choose-2) 
        Delta <- array(0,c(100,n.choose-start+1, ncol(Delta.choose)))
        delta.pwr <- matrix(0,n.choose-start+1, ncol(Delta.choose))
        factor.lev<-as.numeric(unlist(strsplit(input$fl3,",")))
        temp_v<-list()
        temp_c<-list()
        temp_denom<-list()
        temp_n<-list()
        k<-1
        for (n in start:(n.choose)) {
          
          v1=n-1
          
          if (max(input$checkGroup3)==1){
            v <- factor.lev-1
            c <- prod(factor.lev)*n/factor.lev
          } else if (max(input$checkGroup3)==2) {
            v <- (factor.lev-1)%*%t(factor.lev-1)
            v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
            c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
            
          }
          
          
          v.denom <-  prod(factor.lev)*n-1-sum(v)-v1
          
          for (j in 1: length(v) ){
            delta.pwr[(n-n.choose)+nrow(delta.pwr),j]=fsize(input$a3, 1-input$b3, v[j], v.denom, c[j],input$delta_type3);
            
            for (ind in 1: 100){
              
              if(input$a3+1-power[ind]<0.9999)
                (Delta[ind,(n-n.choose)+nrow(delta.pwr),j] <- fsize(input$a3, 1-power[ind], v[j], v.denom, c[j],input$delta_type3))
              else (Delta[ind,(n-n.choose)+nrow(delta.pwr),j]<-NA)
              
            }	}
          temp_n[[k]]<-n
          temp_v[[k]]<-v
          temp_c[[k]]<-c
          temp_denom[[k]]<-v.denom
          k<-k+1
        }
        
        
        
        x<-seq(1:length(full_list1))
        i<-x[full_list1==input$plot_order3] 
        
        
        if(i==max(x))
        {
          plot(Delta[,nrow(delta.pwr),1], power[1:100], type="l", ylab="Power",  
          xlab=ifelse(input$plot_order3=="ALL","Delta",ifelse(input$delta_type3==1,paste0("SD(",input$plot_order3,")/SD(noise)"),paste0("Range(",input$plot_order3,")/SD(noise)"))) ,
          main=ifelse(input$plot_order3=="ALL",paste0("Delta vs Power(=",input$b3,")"),ifelse(input$delta_type3==1,paste0("SD(",input$plot_order3,")/SD(noise) vs Power(=", input$b3,")"),paste0("Range(",input$plot_order3,")/SD(noise) vs Power(=", input$b3,")")))
          ,col=1, lwd=2)
          for(j in 2:(i-1))
          {
            
            points(Delta[,nrow(delta.pwr),j], power[1:100], type="l", col=j,lty=j, lwd=2)
          }
          abline(h=c(0.8,0.9),v=c(1.0,1.5),col="gray")
          
          legend("bottomleft", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
          
          sliderValues <- reactive({})
          
          output$values3 <- renderTable({
            sliderValues()
          })
        }
        
        
        
        else if(i<max(x))
        {     
          plot(Delta[,nrow(delta.pwr),i], power[1:100], type="l", ylab="Power", 
               xlab=ifelse(input$plot_order3=="ALL","Delta",ifelse(input$delta_type3==1,paste0("SD(",input$plot_order3,")/SD(noise)"),paste0("Range(",input$plot_order3,")/SD(noise)"))) ,
               main=ifelse(input$plot_order3=="ALL",paste0("Delta vs Power(=",input$b3,")"),ifelse(input$delta_type3==1,paste0("SD(",input$plot_order3,")/SD(noise) vs Power(=", input$b3,")"),paste0("Range(",input$plot_order3,")/SD(noise) vs Power(=", input$b3,")")))
               ,col=1, lwd=2)
          abline(h=0.8, v= fsize(input$a3,0.2,temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],temp_c[[nrow(delta.pwr)]][i],input$delta_type3) , col="gray", lty=3)
          abline(h=0.9, v =fsize(input$a3,0.1,temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],temp_c[[nrow(delta.pwr)]][i],input$delta_type3) , col="gray", lty=3)
          abline(h= (1-pf(qf((1-input$a3),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]]),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],ncp=ifelse(input$delta_type3==1,(1*(temp_c[[nrow(delta.pwr)]][i]*temp_v[[nrow(delta.pwr)]][i])),
                                                                                                                                                                            (1*temp_c[[nrow(delta.pwr)]][i]/2)))) , v =1.0, col="gray", lty=3)
          abline(h= (1-pf(qf((1-input$a3),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]]),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]],ncp=ifelse(input$delta_type3==1,(1.5^2*(temp_c[[nrow(delta.pwr)]][i]*temp_v[[nrow(delta.pwr)]][i])),
                                                                                                                                                                            (1.5^2*temp_c[[nrow(delta.pwr)]][i]/2)))) , v =1.5 , col="gray", lty=3)
          
          
          if(n.choose>3)
          {
            points(Delta[,nrow(delta.pwr)-1,i], power[1:100], type="l", lty=2, lwd=2)
            
            points(Delta[,nrow(delta.pwr)-2,i], power[1:100], type="l", lty=3, lwd=2)
            
            
          }
          
          else if(n.choose==3)
          {
            points(Delta[,nrow(delta.pwr)-1,i], power[1:100], type="l", lty=2, lwd=2) 
            
            
          }
          
          if(n.choose>3)
          {    
            sliderValues <- reactive({
              
              data.frame(
                R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4),rep(temp_n[[3]],4)),
                Delta = c(round(fsize(input$a3,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type3),3),
                          round(fsize(input$a3,0.1,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type3),3),
                          "1.0","1.5",
                          round(fsize(input$a3,0.2,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type3),3),
                          round(fsize(input$a3,0.1,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type3),3),
                          "1.0","1.5",
                          
                          round(fsize(input$a3,0.2,temp_v[[3]][i],temp_denom[[3]],temp_c[[3]][i],input$delta_type3),3),
                          round(fsize(input$a3,0.1,temp_v[[3]][i],temp_denom[[3]],temp_c[[3]][i],input$delta_type3),3),
                          "1.0","1.5"
                ),
                Power = c("0.8","0.9",
                          round((1-pf(qf((1-input$a3),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type3==1,(1*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                                (1*temp_c[[1]][i]/2)))),3),#(Delta^2)*(c*nu1),
                          
                          round((1-pf(qf((1-input$a3),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type3==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                                (1.5^2*temp_c[[1]][i]/2)))),3),
                          "0.8","0.9",
                          round((1-pf(qf((1-input$a3),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type3==1,(1*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                                (1*temp_c[[2]][i]/2)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a3),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type3==1,(1.5^2*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                                (1.5^2*temp_c[[2]][i]/2)))),3),
                          "0.8","0.9",
                          round((1-pf(qf((1-input$a3),temp_v[[3]][i],temp_denom[[3]]),temp_v[[3]][i],temp_denom[[3]],ncp=ifelse(input$delta_type3==1,(1*(temp_c[[3]][i]*temp_v[[3]][i])),
                                                                                                                                (1*temp_c[[3]][i]/2)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a3),temp_v[[3]][i],temp_denom[[3]]),temp_v[[3]][i],temp_denom[[3]],ncp=ifelse(input$delta_type3==1,(1.5^2*(temp_c[[3]][i]*temp_v[[3]][i])),
                                                                                                                                (1.5^2*temp_c[[3]][i]/2)))),3)),
                stringsAsFactors = FALSE)
              
            })
          }
          else if (n.choose==3)
          {    
            sliderValues <- reactive({
              
              data.frame(
                R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4)),
                Delta = c(round(fsize(input$a3,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type3),3),
                          round(fsize(input$a3,0.1,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type3),3),
                          "1.0","1.5",
                          round(fsize(input$a3,0.2,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type3),3),
                          round(fsize(input$a3,0.1,temp_v[[2]][i],temp_denom[[2]],temp_c[[2]][i],input$delta_type3),3),
                          "1.0","1.5"
                ),
                Power = c("0.8","0.9",
                          round((1-pf(qf((1-input$a3),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type3==1,(1*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                                (1*temp_c[[1]][i]/2)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a3),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type3==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                                (1.5^2*temp_c[[1]][i]/2)))),3),
                          "0.8","0.9",
                          round((1-pf(qf((1-input$a3),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type3==1,(1*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                                (1*temp_c[[2]][i]/2)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a3),temp_v[[2]][i],temp_denom[[2]]),temp_v[[2]][i],temp_denom[[2]],ncp=ifelse(input$delta_type3==1,(1.5^2*(temp_c[[2]][i]*temp_v[[2]][i])),
                                                                                                                                (1.5^2*temp_c[[2]][i]/2)))),3)
                          
                ),
                stringsAsFactors = FALSE)
              
            })
          }
          else if (n.choose==2)
          {    
            sliderValues <- reactive({
              
              data.frame(
                R=c(rep(temp_n[[1]],4) ),
                Delta = c(round(fsize(input$a3,0.2,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type3),3),
                          round(fsize(input$a3,0.1,temp_v[[1]][i],temp_denom[[1]],temp_c[[1]][i],input$delta_type3),3),
                          "1.0","1.5" 
                ),
                Power = c("0.8","0.9",
                          round((1-pf(qf((1-input$a3),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type3==1,(1*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                                (1*temp_c[[1]][i]/2)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a3),temp_v[[1]][i],temp_denom[[1]]),temp_v[[1]][i],temp_denom[[1]],ncp=ifelse(input$delta_type3==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i])),
                                                                                                                                (1.5^2*temp_c[[1]][i]/2)))),3) 
                          
                ),
                stringsAsFactors = FALSE)
              
            })
          }
          # Show the values in an HTML table ----
          output$values3 <- renderTable({
            sliderValues()
          })
        }
        
      })
      
      
      
      
      output$Size_graph3<-renderPlot({
        FF2<-sampleSize.RCBD(input$nf3, as.numeric(unlist(strsplit(input$fl3,","))),  order=max(input$checkGroup3), delta_type=input$delta_type3,Deltao=ifelse(rep(input$delta_type3==1,3),c(input$de1_3,input$de2_3,input$de3_3),c(input$de11_3,input$de12_3,input$de13_3)), beta=1-input$b3, alpha=input$a3)
        
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        
        factor.lev<-as.numeric(unlist(strsplit(input$fl3,",")))
        
        Delta <- matrix(0,100,ncol(Delta.choose)) 
        
        for (n in 2:100) {
          
          v1=n-1
          
          if (max(input$checkGroup3)==1){
            v <- factor.lev-1
            c <- prod(factor.lev)*n/factor.lev
          } else if (max(input$checkGroup3)==2) {
            v <- (factor.lev-1)%*%t(factor.lev-1)
            v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
            c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
            
          }
          
          
          v.denom <-  prod(factor.lev)*n-1-sum(v)-v1
          
          for (j in 1: length(v) ){
            Delta[n,j] <- fsize(input$a3, 1-input$b3, v[j], v.denom, c[j],input$delta_type3)
          }
        }
        
        
        plot(2:100, Delta[2:100,1], type="l", xlim=c(0,min(100,n.choose+5)), ylim=c(0,max(ifelse(rep(input$delta_type3==1,2),c(input$de1_3/input$de3_3,input$de2_3/input$de3_3),c(input$de11_3/input$de13_3,input$de12_3/input$de13_3)))*1.5), ylab="Delta", xlab="The number of replications", 
             main="The number of replications(r) vs Delta",col=1, lwd=2)
        for (i in 2:ncol(Delta)) 
          lines(2:100, Delta[2:100,i], type="l", lty=i, lwd=2,col=i)
        
        abline(h=max(ifelse(rep(input$delta_type3==1,2),c(input$de1_3/input$de3_3,input$de2_3/input$de3_3),c(input$de11_3/input$de13_3,input$de12_3/input$de13_3))), v=FF2$n,col="gray", lty=3)
        legend("top", legend=paste0("power=", input$b3), adj=0, bty="n")
        
        legend("topright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
        
        
      })
      
      
      output$power_graph3<-renderPlot({
        FF2<-sampleSize.RCBD(input$nf3, as.numeric(unlist(strsplit(input$fl3,","))),  order=max(input$checkGroup3), delta_type=input$delta_type3,Deltao=ifelse(rep(input$delta_type3==1,3),c(input$de1_3,input$de2_3,input$de3_3),c(input$de11_3,input$de12_3,input$de13_3)), beta=1-input$b3, alpha=input$a3)
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        Deltao <- c(1,1.5,2)
        Delta <- pwr <- array(0,c(100, ncol(Delta.choose),3))
        delta.pwr <- array(0,c(100,ncol(Delta.choose),3))
        factor.lev<-as.numeric(unlist(strsplit(input$fl3,",")))
        
        
        for (deltao in 1: 3){
          for (n in 2:100) {
            
            v1=n-1
            
            if (max(input$checkGroup3)==1){
              v <- factor.lev-1
              c <- prod(factor.lev)*n/factor.lev
            } else if (max(input$checkGroup3)==2) {
              v <- (factor.lev-1)%*%t(factor.lev-1)
              v <- c(factor.lev-1,v[upper.tri(v, diag=FALSE)])
              c <- prod(factor.lev)*n/c(factor.lev, (factor.lev%*%t(factor.lev))[upper.tri((factor.lev)%*%t(factor.lev), diag=FALSE)])
              
            }
            
            
            v.denom <- prod(factor.lev)*n-1-sum(v)-v1
            
            for (j in 1: length(v)){ 
              pwr[n,j,deltao]<-round((1-pf(qf((1-input$a3),v[j],v.denom),v[j],v.denom,ncp=ifelse(input$delta_type3==1,(Deltao[deltao]^2*(c[j]*v[j])),
                                                                                                 (Deltao[deltao]^2*c[j]/2)))),3)
              
              
            } 
            
            
            
            
          }
        }
        
        x<-c(1,2,3)
        j<-x[Deltao==input$plot_delta3]
        
        plot(2:100, pwr[2:100,1, j], type="l", ylim=c(0,1), xlim=c(0,min(100,n.choose+5)), ylab="Power", xlab="The number of replications", 
             main="The number of replications(r) vs Power",col=1, lwd=2)
        for (i in 2:ncol(Delta)) 
          lines(2:100, pwr[2:100,i,j], type="l", lty=i, lwd=2,col=i)
        
        abline(h=0.8, v=FF2$n,col="gray", lty=3)
        
        abline(h=0.9, col="gray", lty=3)
        legend("bottomright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
        
      }
      
      
      )
    })
    
    
    
    observeEvent(input$do4, {
      #fsize.m
      w_main_list<-list()
      for(i in 1:input$wf)
      {
        w_main_list[[i]]<-toupper(letters)[i]
      }
      s_main_list<-list()
      for(i in 1:input$sf)
      {
        s_main_list[[i]]<-toupper(letters)[length(w_main_list)+i]
      }
      
      if(max(input$checkGroup4)==1)
      {
        
        full_list1<-c(w_main_list,s_main_list,"ALL")
        full_list<-c(w_main_list,s_main_list)
        updateSelectInput(session,"plot_order4",choices=full_list1)
      }
      else if(max(input$checkGroup4)==2)
      {
        
        if(input$sf>1)
        {
        k<-1
        s_two_list<-list()
        
        for(i in 1:(length(s_main_list)-1)) {
          for(j in 1:(length(s_main_list)-i)){
            s_two_list[[k]]<-paste(s_main_list[[i]],"*",s_main_list[[i+j]])
            k<-k+1
          }
        }
        
      }
      
        if(input$wf>1)
        {
        k<-1
        w_two_list<-list()
        
        for(i in 1:(length(w_main_list)-1)) {
          for(j in 1:(length(w_main_list)-i)){
            w_two_list[[k]]<-paste(w_main_list[[i]],"*",w_main_list[[i+j]])
            k<-k+1
          }
        }
        }
        k<-1
        ws_two_list<-list()
        
        for(i in 1:length(w_main_list)){
          for (j in 1:length(s_main_list)){
            ws_two_list[[k]]<-paste(w_main_list[[i]],"*",s_main_list[[j]])
            k<-k+1
          }
        }
        
        if(input$sf>1 & input$wf>1)
        {
        full_list1<-c(w_main_list,w_two_list,s_main_list,s_two_list,ws_two_list,"ALL")
        full_list<-c(w_main_list,w_two_list,s_main_list,s_two_list,ws_two_list)
        }
        
        else if(input$sf>1) {
        full_list1<-c(w_main_list,s_main_list,s_two_list,ws_two_list,"ALL")
        full_list<-c(w_main_list,s_main_list,s_two_list,ws_two_list)
        }
        
        else if(input$wf>1) {
          full_list1<-c(w_main_list,w_two_list,s_main_list, ws_two_list,"ALL")
          full_list<-c(w_main_list,w_two_list,s_main_list, ws_two_list)
        }
        else {
          
          full_list1<-c(w_main_list, s_main_list, ws_two_list,"ALL")
          full_list<-c(w_main_list, s_main_list, ws_two_list)
        }
        
        updateSelectInput(session,"plot_order4",choices= full_list1)
      }
      
      fsize <-  function(alpha, beta, nu1, nu2, c,delta_type){
        fc <- qf(1-alpha, nu1, nu2)
        fl <- ncparamF(alpha, beta, nu1, nu2)/2
        if (delta_type==1) (Delta<-sqrt(2*fl/(c*nu1))) #ncp<-Delta^2*c*nu1
        else if (delta_type==2) (Delta<-sqrt(4*fl/c)) #ncp<-Delta^2*c/2
        
        return(Delta)
      }
      
      list_fc<-function(list)
      {
        x<- length(list) 
        
        list1<-list[[1]]
        for( i in 1:(x-1))
          list1<-paste(list1,list[[i+1]],sep="+")
        return(list1=list1)
      }
      
      sampleSize.split<- function(whole.factor, whole.factor.lev, split.factor, split.factor.lev, delta_type=1,order=c(1,2), Deltao=c(1,1,1,1), alpha=0.05, beta=0.2){
        main_n<-0
        two_n<-0
        nn<-0
        whole.Delta <- NULL
        splot.Delta<-NULL
        
        
        for (n in 2:100){
          if (order==1){
            
            v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
            c.whole <- NULL ; c.split <- NULL
            whole.Delta <- NULL; split.Delta <- NULL
            
            v.rep <- n-1
            
            #whole
            v.whole <- whole.factor.lev-1
            v.whole.denom <- prod(whole.factor.lev)*n-1-sum(whole.factor.lev-1)-v.rep
            
            c.whole <- prod(whole.factor.lev)*prod(split.factor.lev)*n/whole.factor.lev
            
            for (i in 1: length(v.whole)){
              whole.Delta[i] <- fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i],delta_type)*sqrt(prod(split.factor.lev)+1)
            }	
            
            #split
            v.split <- split.factor.lev-1
            v.split.denom <- prod(whole.factor.lev)*prod(split.factor.lev)*n-1-(prod(whole.factor.lev)*n-1)-sum(split.factor.lev-1)
            
            c.split <- prod(whole.factor.lev)*prod(split.factor.lev)*n/split.factor.lev
            
            for (i in 1: length(v.split)){
              split.Delta[i] <- fsize(alpha, beta, v.split[i], v.split.denom, c.split[i],delta_type)
            }
            
            if ( max(whole.Delta)<=Deltao[1]/Deltao[3] & max(split.Delta)<=Deltao[1]/Deltao[4]  ) 
              (nn<-n)
            
          }
          
          else if (order==2) {
            
            
            
            v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
            c.whole <- NULL ; c.split <- NULL
            whole.Delta <- NULL; split.Delta <- NULL
            
            v.rep <- n-1
            
            #whole
            v.whole <- (whole.factor.lev-1)%*%t(whole.factor.lev-1)
            v.whole <- c(whole.factor.lev-1,v.whole[upper.tri(v.whole, diag=FALSE)])
            v.whole.denom <- prod(whole.factor.lev)*n-1-sum(v.whole)-v.rep
            
            c.whole <- prod(whole.factor.lev)*prod(split.factor.lev)*n/c(whole.factor.lev, (whole.factor.lev%*%t(whole.factor.lev))[upper.tri((whole.factor.lev)%*%t(whole.factor.lev), diag=FALSE)])
            
            
            
            #split
            v.split <- (split.factor.lev-1)%*%t(split.factor.lev-1)
            v.split.temp <- (whole.factor.lev-1)%*%t(split.factor.lev-1)
            v.split <- c(split.factor.lev-1, v.split[upper.tri(v.split, diag=FALSE)], as.vector(t(v.split.temp)))
            
            v.split.denom <- prod(whole.factor.lev)*prod(split.factor.lev)*n-prod(whole.factor.lev)*n-sum(v.split)
            
            c.split <- prod(whole.factor.lev)*prod(split.factor.lev)*n/c(split.factor.lev, (split.factor.lev%*%t(split.factor.lev))[upper.tri((split.factor.lev)%*%t(split.factor.lev), diag=FALSE)], as.vector(t(whole.factor.lev%*%t(split.factor.lev))))
            
            if(main_n==0)
            {
              for (i in 1:whole.factor){
                whole.Delta[i] <-  fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i],delta_type)*sqrt(prod(split.factor.lev)+1)
              }
              for (i in 1:split.factor){
                split.Delta[i] <-fsize(alpha, beta, v.split[i], v.split.denom, c.split[i],delta_type)
              }
              if (max(whole.Delta[1:whole.factor])<=Deltao[1]/Deltao[3] & max(split.Delta[1:split.factor])<=Deltao[1]/Deltao[4] ) (main_n<-n)
            }
            
            if(two_n==0)
            {
              for(i in (split.factor+1) : length(v.split)){
                split.Delta[i] <- fsize(alpha, beta, v.split[i], v.split.denom, c.split[i],delta_type)
              }
              
              if(whole.factor>1){
                for (i in (whole.factor+1): length(v.whole)){
                  whole.Delta[i] <- fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i],delta_type)*sqrt(prod(split.factor.lev)+1)
                }	
                
                if (max(whole.Delta[(whole.factor+1): length(v.whole)])<=Deltao[2]/Deltao[3] & max(split.Delta[(split.factor+1) : length(v.split)])<=Deltao[2]/Deltao[4] ) (two_n<-n)
              }
              else 
                ( if (  max(split.Delta[(split.factor+1) : length(v.split)])<=Deltao[2]/Deltao[4] ) (two_n<-n))
            }
            if(main_n >0 & two_n>0) (nn<-max(main_n,two_n))
          }
          
          
          
          if(nn>0) 
          {
            for(i in 1:length(v.whole)){
              whole.Delta[i]<-fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i],delta_type)*sqrt(prod(split.factor.lev)+1)
            }
            for(i in 1:length(v.split)){
              split.Delta[i] <- fsize(alpha, beta, v.split[i], v.split.denom, c.split[i],delta_type)
            }
            
            Delta <- c(whole.Delta, split.Delta)
            break
          }
        }
        
        return(list(n=nn, Delta=Delta))	
      }
           
      output$Size1_4<-renderText({sampleSize.split(input$wf, as.numeric(unlist(strsplit(input$wfl,","))), input$sf, as.numeric(unlist(strsplit(input$sfl,","))),  
                                                   delta_type=input$delta_type4,order=max(input$checkGroup4), Deltao=ifelse(rep(input$delta_type4==1,4),c(input$de1_4,input$de2_4,input$de3_4,input$de4_4),c(input$de11_4,input$de12_4,input$de13_4,input$de14_4)), beta=1-input$b4, alpha=input$a4)[[1]]}) 
      output$Size2_4<-renderText({sampleSize.split(input$wf, as.numeric(unlist(strsplit(input$wfl,","))), input$sf, as.numeric(unlist(strsplit(input$sfl,","))),  
                                                   delta_type=input$delta_type4,order=max(input$checkGroup4), Deltao=ifelse(rep(input$delta_type4==1,4),c(input$de1_4,input$de2_4,input$de3_4,input$de4_4),c(input$de11_4,input$de12_4,input$de13_4,input$de14_4)), beta=1-input$b4, alpha=input$a4)[[2]]}) 
      
      output$list1_4<-renderText({list_fc(full_list)})
      
      output$Delta_graph4<-renderPlot({
        FF2<-sampleSize.split(input$wf, as.numeric(unlist(strsplit(input$wfl,","))), input$sf, as.numeric(unlist(strsplit(input$sfl,","))),  
                              delta_type=input$delta_type4,order=max(input$checkGroup4), Deltao=ifelse(rep(input$delta_type4==1,4),c(input$de1_4,input$de2_4,input$de3_4,input$de4_4),c(input$de11_4,input$de12_4,input$de13_4,input$de14_4)), beta=1-input$b4, alpha=input$a4)
        
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        
        start <- ifelse((n.choose-1)<=2, 2, n.choose-2) 
        Delta <- array(0,c(100,n.choose-start+1, ncol(Delta.choose)))
        delta.pwr <- matrix(0,n.choose-start+1, ncol(Delta.choose)) 
        wfl<-as.numeric(unlist(strsplit(input$wfl,",")))
        sfl<-as.numeric(unlist(strsplit(input$sfl,",")))
        temp_n<-list()
        temp_v<-list()
        temp_c<-list()
        temp_denom<-list() 
        
        k<-1
        for (n in start:(n.choose)) {
          
          if (max(input$checkGroup4)==1){
            
            v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
            c.whole <- NULL ; c.split <- NULL
            whole.Delta <- NULL; split.Delta <- NULL
            
            v.rep <- n-1
            
            #whole
            v.whole <-  wfl-1
            v.whole.denom <- prod( wfl)*n-1-sum( wfl-1)-v.rep
            
            c.whole <- prod( wfl)*prod( sfl)*n/ wfl
            
            for (i in 1: length(v.whole)){
              whole.Delta[i] <- fsize(input$a4, 1-input$b4, v.whole[i], v.whole.denom, c.whole[i],input$delta_type4)*sqrt(prod( sfl)+1)
              
              for (ind in 1: 100){
                if(input$a4+1-power[ind]<0.9999)
                  (Delta[ind,(n-n.choose)+nrow(delta.pwr),i] <- fsize(input$a, 1-power[ind],v.whole[i], v.whole.denom, c.whole[i],input$delta_type4)*sqrt(prod( sfl)+1))
                else ( Delta[ind,(n-n.choose)+nrow(delta.pwr),i]<-NA)
              }
            }	
            
            #split
            v.split <-  sfl-1
            v.split.denom <- prod( wfl)*prod( sfl)*n-1-(prod( wfl)*n-1)-sum( sfl-1)
            
            c.split <- prod( wfl)*prod( sfl)*n/ sfl
            
            for (i in 1: length(v.split)){
              split.Delta[i] <- fsize(input$a4, 1-input$b4, v.split[i], v.split.denom, c.split[i],input$delta_type4)
              
              for (ind in 1: 100){
                if(input$a4+1-power[ind]<0.9999)
                  (Delta[ind,(n-n.choose)+nrow(delta.pwr),(length(v.whole)+i)] <- fsize(input$a4, 1-power[ind], v.split[i], v.split.denom, c.split[i],input$delta_type4))
                else ( Delta[ind,(n-n.choose)+nrow(delta.pwr),(length(v.whole)+i)] <-NA)
              }
            }
            
            
            delta.pwr[(n-n.choose)+nrow(delta.pwr),1:(length(v.whole))] <- whole.Delta
            delta.pwr[(n-n.choose)+nrow(delta.pwr),(1+length(v.whole)):(length(v.whole)+length(v.split))] <- split.Delta
          } 
          else if (max(input$checkGroup4)==2){
            
            v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
            c.whole <- NULL ; c.split <- NULL
            whole.Delta <- NULL; split.Delta <- NULL
            
            v.rep <- n-1
            
            #whole
            v.whole <- ( wfl-1)%*%t( wfl-1)
            v.whole <- c( wfl-1,v.whole[upper.tri(v.whole, diag=FALSE)])
            v.whole.denom <- prod( wfl)*n-1-sum(v.whole)-v.rep
            
            c.whole <- prod( wfl)*prod( sfl)*n/c( wfl, ( wfl%*%t( wfl))[upper.tri(( wfl)%*%t( wfl), diag=FALSE)])
            
            for (i in 1: length(v.whole)){
              whole.Delta[i] <- fsize(input$a4, 1-input$b4, v.whole[i], v.whole.denom, c.whole[i],input$delta_type4)*sqrt(prod( sfl)+1)
              
              for (ind in 1: 100){
                
                if(input$a4+1-power[ind]<0.9999)
                  (Delta[ind,(n-n.choose)+nrow(delta.pwr),i] <- fsize(input$a4, 1-power[ind],v.whole[i], v.whole.denom, c.whole[i],input$delta_type4)*sqrt(prod( sfl)+1))
                else (Delta[ind,(n-n.choose)+nrow(delta.pwr),i] <-NA)
              }
            }	
            
            #split
            v.split <- ( sfl-1)%*%t( sfl-1)
            v.split.temp <- (wfl-1)%*%t( sfl-1)
            v.split <- c( sfl-1, v.split[upper.tri(v.split, diag=FALSE)], as.vector(t(v.split.temp)))
            
            v.split.denom <- prod( wfl)*prod(sfl)*n-prod( wfl)*n-sum(v.split)
            
            c.split <- prod( wfl)*prod( sfl)*n/c( sfl, ( sfl%*%t( sfl))[upper.tri(( sfl)%*%t(sfl), diag=FALSE)], as.vector(t( wfl%*%t( sfl))))
            
            for (i in 1: length(v.split)){
              split.Delta[i] <- fsize(input$a4, 1-input$b4, v.split[i], v.split.denom, c.split[i],input$delta_type4)
              
              for (ind in 1: 100){
                
                if(input$a4+1-power[ind]<0.9999)
                  (Delta[ind,(n-n.choose)+nrow(delta.pwr),(length(v.whole)+i)] <- fsize(input$a4, 1-power[ind], v.split[i], v.split.denom, c.split[i],input$delta_type4))
                else (Delta[ind,(n-n.choose)+nrow(delta.pwr),(length(v.whole)+i)] <-NA)
              }
            }
            
            
            delta.pwr[(n-n.choose)+nrow(delta.pwr),1:length(v.whole)] <- whole.Delta
            delta.pwr[(n-n.choose)+nrow(delta.pwr),(1+length(v.whole)):(length(v.whole)+length(v.split))] <- split.Delta;
          }
          
          temp_n[[k]]<-n
          temp_v[[k]]<-c(v.whole,v.split) 
          temp_c[[k]]<-c(c.whole,c.split)
          temp_denom[[k]]<-c(v.whole.denom,v.split.denom) 
          
          k<-k+1
        }
        
        
        
        
        
        x<-seq(1:length(full_list1))
        i<-x[full_list1==input$plot_order4] 
        
        if(i==max(x))
        {
          plot(Delta[,nrow(delta.pwr),1], power[1:100], type="l", ylab="Power",  
          xlab=ifelse(input$plot_order4=="ALL","Delta",ifelse(input$delta_type4==1,paste0("SD(",input$plot_order4,")/SD(noise)"),paste0("Range(",input$plot_order4,")/SD(noise)"))) ,
          main=ifelse(input$plot_order4=="ALL",paste0("Delta vs Power(=",input$b4,")"),ifelse(input$delta_type4==1,paste0("SD(",input$plot_order4,")/SD(noise) vs Power(=", input$b4,")"),paste0("Range(",input$plot_order4,")/SD(noise) vs Power(=", input$b4,")")))
          ,col=1, lwd=2)
          for(j in 2:(i-1))
          {
            
            points(Delta[,nrow(delta.pwr),j], power[1:100], type="l", col=j,lty=j, lwd=2)
          }
          abline(h=c(0.8,0.9),v=c(1.0,1.5),col="gray")
          
          legend("bottomleft", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
          
          sliderValues <- reactive({})
          
          output$values4 <- renderTable({
            sliderValues()
          })
        }
        
        
        
        else if(i<max(x))
        {     
          plot(Delta[,nrow(delta.pwr),i], power[1:100], type="l", ylab="Power", xlab=ifelse(input$plot_order4=="ALL","Delta",ifelse(input$delta_type4==1,paste0("SD(",input$plot_order4,")/SD(noise)"),paste0("Range(",input$plot_order4,")/SD(noise)"))) ,
               main=ifelse(input$plot_order4=="ALL",paste0("Delta vs Power(=",input$b4,")"),ifelse(input$delta_type4==1,paste0("SD(",input$plot_order4,")/SD(noise) vs Power(=", input$b4,")"),paste0("Range(",input$plot_order4,")/SD(noise) vs Power(=", input$b4,")")))
               ,col=1, lwd=2)
          abline(h=0.8, v= fsize(input$a4,0.2,temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]][ifelse(i<=length(v.whole),1,2)],temp_c[[nrow(delta.pwr)]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1) , col="gray", lty=3)
          abline(h=0.9, v =fsize(input$a4,0.1,temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]][ifelse(i<=length(v.whole),1,2)],temp_c[[nrow(delta.pwr)]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1) , col="gray", lty=3)
          abline(h= (1-pf(qf((1-input$a4),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]][ifelse(i<=length(v.whole),1,2)]),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1*(temp_c[[nrow(delta.pwr)]][i]*temp_v[[nrow(delta.pwr)]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                                                            (1*temp_c[[nrow(delta.pwr)]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))) , v =1.0, col="gray", lty=3)
          abline(h= (1-pf(qf((1-input$a),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]][ifelse(i<=length(v.whole),1,2)]),temp_v[[nrow(delta.pwr)]][i],temp_denom[[nrow(delta.pwr)]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1.5^2*(temp_c[[nrow(delta.pwr)]][i]*temp_v[[nrow(delta.pwr)]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                                                           (1.5^2*temp_c[[nrow(delta.pwr)]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))) , v =1.5 , col="gray", lty=3)
          
          
          if(n.choose>3)
          {
            points(Delta[,nrow(delta.pwr)-1,i], power[1:100], type="l", lty=2, lwd=2)
            
            points(Delta[,nrow(delta.pwr)-2,i], power[1:100], type="l", lty=3, lwd=2)
            
            
          }
          
          else if(n.choose==3)
          {
            points(Delta[,nrow(delta.pwr)-1,i], power[1:100], type="l", lty=2, lwd=2) 
            
            
          }
          
          if(n.choose>3)
          {    
            sliderValues <- reactive({
              
              data.frame(
                R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4),rep(temp_n[[3]],4)),
                Delta = c(round(fsize(input$a4,0.2,temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],temp_c[[1]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          round(fsize(input$a4,0.1,temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],temp_c[[1]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          "1.0","1.5",
                          round(fsize(input$a4,0.2,temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],temp_c[[2]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          round(fsize(input$a4,0.1,temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],temp_c[[2]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          "1.0","1.5",
                          
                          round(fsize(input$a4,0.2,temp_v[[3]][i],temp_denom[[3]][ifelse(i<=length(v.whole),1,2)],temp_c[[3]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          round(fsize(input$a4,0.1,temp_v[[3]][i],temp_denom[[3]][ifelse(i<=length(v.whole),1,2)],temp_c[[3]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          "1.0","1.5"
                ),
                Power = c("0.8","0.9",
                          round((1-pf(qf((1-input$a4),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)]),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1*(temp_c[[1]][i]*temp_v[[1]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1*temp_c[[1]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),#(Delta^2)*(c*nu1),
                          
                          round((1-pf(qf((1-input$a4),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)]),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1.5^2*temp_c[[1]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),
                          "0.8","0.9",
                          round((1-pf(qf((1-input$a4),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)]),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1*(temp_c[[2]][i]*temp_v[[2]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1*temp_c[[2]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a4),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)]),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1.5^2*(temp_c[[2]][i]*temp_v[[2]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1.5^2*temp_c[[2]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),
                          "0.8","0.9",
                          round((1-pf(qf((1-input$a4),temp_v[[3]][i],temp_denom[[3]][ifelse(i<=length(v.whole),1,2)]),temp_v[[3]][i],temp_denom[[3]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1*(temp_c[[3]][i]*temp_v[[3]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1*temp_c[[3]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a4),temp_v[[3]][i],temp_denom[[3]][ifelse(i<=length(v.whole),1,2)]),temp_v[[3]][i],temp_denom[[3]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1.5^2*(temp_c[[3]][i]*temp_v[[3]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1.5^2*temp_c[[3]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3)),
                stringsAsFactors = FALSE)
              
            })
          }
          else if (n.choose==3)
          {    
            sliderValues <- reactive({
              
              data.frame(
                R=c(rep(temp_n[[1]],4),rep(temp_n[[2]],4)),
                Delta = c(round(fsize(input$a4,0.2,temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],temp_c[[1]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          round(fsize(input$a4,0.1,temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],temp_c[[1]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          "1.0","1.5",
                          round(fsize(input$a4,0.2,temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],temp_c[[2]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          round(fsize(input$a4,0.1,temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],temp_c[[2]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          "1.0","1.5"
                ),
                Power = c("0.8","0.9",
                          round((1-pf(qf((1-input$a4),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)]),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1*(temp_c[[1]][i]*temp_v[[1]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1*temp_c[[1]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a4),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)]),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1.5^2*temp_c[[1]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),
                          "0.8","0.9",
                          round((1-pf(qf((1-input$a4),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)]),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1*(temp_c[[2]][i]*temp_v[[2]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1*temp_c[[2]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a4),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)]),temp_v[[2]][i],temp_denom[[2]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1.5^2*(temp_c[[2]][i]*temp_v[[2]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                (1.5^2*temp_c[[2]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3)
                          
                ),
                stringsAsFactors = FALSE)
              
            })
          }
          else if (n.choose==2)
          {    
            sliderValues <- reactive({
              
              data.frame(
                R=c(rep(temp_n[[1]],4) ),
                Delta = c(round(fsize(input$a4,0.2,temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],temp_c[[1]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          round(fsize(input$a4,0.1,temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],temp_c[[1]][i],input$delta_type4)*ifelse(i<=length(v.whole),sqrt(prod( sfl)+1),1),3),
                          "1.0","1.5" 
                ),
                Power = c("0.8","0.9",
                          round((1-pf(qf((1-input$a4),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)]),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1*(temp_c[[1]][i]*temp_v[[1]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                                                                                (1*temp_c[[1]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3),#(Delta^2)*(c*nu1)
                          round((1-pf(qf((1-input$a4),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)]),temp_v[[1]][i],temp_denom[[1]][ifelse(i<=length(v.whole),1,2)],ncp=ifelse(input$delta_type4==1,(1.5^2*(temp_c[[1]][i]*temp_v[[1]][i]))/ifelse(i<=length(v.whole),prod(sfl)+1,1),
                                                                                                                                                                                                (1.5^2*temp_c[[1]][i]/2)/ifelse(i<=length(v.whole),prod(sfl)+1,1)))),3) 
                          
                ),
                stringsAsFactors = FALSE)
              
            })
          }
          # Show the values in an HTML table ----
          output$values4 <- renderTable({
            sliderValues()
          })
        }
        
      })
          
          
      
      output$Size_graph4<-renderPlot({
        FF2<-sampleSize.split(input$wf, as.numeric(unlist(strsplit(input$wfl,","))), input$sf, as.numeric(unlist(strsplit(input$sfl,","))),  
                              delta_type=input$delta_type4,order=max(input$checkGroup4), Deltao=ifelse(rep(input$delta_type4==1,4),c(input$de1_4,input$de2_4,input$de3_4,input$de4_4),c(input$de11_4,input$de12_4,input$de13_4,input$de14_4)), beta=1-input$b4, alpha=input$a4)
        
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        wfl<-as.numeric(unlist(strsplit(input$wfl,",")))
        sfl<-as.numeric(unlist(strsplit(input$sfl,",")))
        Delta <- matrix(0,100,ncol(Delta.choose)) 
        
        for (n in 2:100) {
          
          if (max(input$checkGroup4)==1){
            
            v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
            c.whole <- NULL ; c.split <- NULL 
            
            v.rep <- n-1
            
            #whole
            v.whole <-  wfl-1
            v.whole.denom <- prod( wfl)*n-1-sum( wfl-1)-v.rep
            
            c.whole <- prod( wfl)*prod( sfl)*n/ wfl
            
            for (i in 1: length(v.whole)){
              Delta[n,i] <- fsize(input$a4, 1-input$b4, v.whole[i], v.whole.denom, c.whole[i],input$delta_type4)*sqrt(prod( sfl)+1)
            }	
            
            #split
            v.split <-  sfl-1
            v.split.denom <- prod( wfl)*prod( sfl)*n-1-(prod( wfl)*n-1)-sum( sfl-1)
            
            c.split <- prod( wfl)*prod( sfl)*n/ sfl
            
            for (i in 1: length(v.split)){
              Delta[n,(length(v.whole)+i)] <- fsize(input$a4, 1-input$b4, v.split[i], v.split.denom, c.split[i],input$delta_type4)
              
            }
            
            
          } 
          else if (max(input$checkGroup4)==2){
            
            v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
            c.whole <- NULL ; c.split <- NULL 
            
            v.rep <- n-1
            
            #whole
            v.whole <- ( wfl-1)%*%t( wfl-1)
            v.whole <- c( wfl-1,v.whole[upper.tri(v.whole, diag=FALSE)])
            v.whole.denom <- prod( wfl)*n-1-sum(v.whole)-v.rep
            
            c.whole <- prod( wfl)*prod( sfl)*n/c( wfl, ( wfl%*%t( wfl))[upper.tri(( wfl)%*%t( wfl), diag=FALSE)])
            
            for (i in 1: length(v.whole)){
              Delta[n,i] <- fsize(input$a4, 1-input$b4, v.whole[i], v.whole.denom, c.whole[i],input$delta_type4)*sqrt(prod( sfl)+1)
              
            }	
            
            #split
            v.split <- ( sfl-1)%*%t( sfl-1)
            v.split.temp <- (wfl-1)%*%t( sfl-1)
            v.split <- c( sfl-1, v.split[upper.tri(v.split, diag=FALSE)], as.vector(t(v.split.temp)))
            
            v.split.denom <- prod( wfl)*prod(sfl)*n-prod( wfl)*n-sum(v.split)
            
            c.split <- prod( wfl)*prod( sfl)*n/c( sfl, ( sfl%*%t( sfl))[upper.tri(( sfl)%*%t(sfl), diag=FALSE)], as.vector(t( wfl%*%t( sfl))))
            
            for (i in 1: length(v.split)){
              Delta[n,(length(v.whole)+i)] <- fsize(input$a4, 1-input$b4, v.split[i], v.split.denom, c.split[i],input$delta_type4)
              
            }
            
          }
          
          
        }
        
        
        plot(2:100, Delta[2:100,1], type="l", xlim=c(0,min(100,n.choose+5)), ylim=c(0,max(ifelse(rep(input$delta_type4==1,2),c(input$de1_4/min(input$de3_4,input$de4_4),input$de2_4/min(input$de3_4,input$de4_4)),c(input$de11_4/min(input$de13_4,input$de14_4),input$de12_4/min(input$de13_4,input$de14_4))))*1.5), ylab="Delta", xlab="The number of replications", 
             main="The number of replications(r) vs Delta",col=1, lwd=2)
        for (i in 2:ncol(Delta)) 
          lines(2:100, Delta[2:100,i], type="l", lty=i, lwd=2,col=i)
        
        abline(h=max(ifelse(rep(input$delta_type4==1,2),c(input$de1_4/min(input$de3_4,input$de4_4),input$de2_4/min(input$de3_4,input$de4_4)),c(input$de11_4/min(input$de13_4,input$de14_4),input$de12_4/min(input$de13_4,input$de14_4)))), v=FF2$n,col="gray", lty=3)
        legend("top", legend=paste0("power=", input$b4), adj=0, bty="n")
        
        legend("topright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
        
        
      })
      #################################
      
      
      output$power_graph4<-renderPlot({
        FF2<-sampleSize.split(input$wf, as.numeric(unlist(strsplit(input$wfl,","))), input$sf, as.numeric(unlist(strsplit(input$sfl,","))),  
                              delta_type=input$delta_type4,order=max(input$checkGroup4), Deltao=ifelse(rep(input$delta_type4==1,4),c(input$de1_4,input$de2_4,input$de3_4,input$de4_4),c(input$de11_4,input$de12_4,input$de13_4,input$de14_4)), beta=1-input$b4, alpha=input$a4)
        (n.choose <- FF2$n); 
        (Delta.choose <- data.frame(t(FF2$Delta)))
        gsize<-length(Delta.choose)
        power <- round(seq(0,1,length.out=101),3)
        Deltao <- c(1,1.5,2)
        Delta <- pwr <- array(0,c(100, ncol(Delta.choose),3))
        delta.pwr <- array(0,c(100,ncol(Delta.choose),3))
        wfl<-as.numeric(unlist(strsplit(input$wfl,",")))
        sfl<-as.numeric(unlist(strsplit(input$sfl,",")))
        
        
        for (deltao in 1: 3){
          
          for (n in 2:100) {
            
            if (max(input$checkGroup4)==1){
              
              v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
              c.whole <- NULL ; c.split <- NULL 
              
              v.rep <- n-1
              
              #whole
              v.whole <-  wfl-1
              v.whole.denom <- prod( wfl)*n-1-sum( wfl-1)-v.rep
              
              c.whole <- prod( wfl)*prod( sfl)*n/ wfl
              
              for (i in 1: length(v.whole)){
                pwr[n,i,deltao]<- (1-pf(qf((1-input$a4),v.whole[i], v.whole.denom),v.whole[i],v.whole.denom,ncp=ifelse(input$delta_type4==1,(Deltao[deltao]^2*(c.whole[i]*v.whole[i]))/(prod(sfl)+1),
                                                                                                                        (Deltao[deltao]^2*c.whole[i]/2)/(prod(sfl)+1))))#(Delta^2)*(c*nu1)
              }	
              
              #split
              v.split <-  sfl-1
              v.split.denom <- prod( wfl)*prod( sfl)*n-1-(prod( wfl)*n-1)-sum( sfl-1)
              
              c.split <- prod( wfl)*prod( sfl)*n/ sfl
              
              for (i in 1: length(v.split)){
                pwr[n,(i+length(v.whole)),deltao]<- (1-pf(qf((1-input$a4),v.split[i], v.split.denom),v.split[i],v.split.denom,ncp=ifelse(input$delta_type4==1,(Deltao[deltao]^2*(c.split[i]*v.split[i])) ,
                                                                                                                                          (Deltao[deltao]^2*c.split[i]/2) )))#(Delta^2)*(c*nu1)
              }
              
              
              
            } 
            else if (max(input$checkGroup4)==2){
              
              v.whole <- NULL ; v.split <- NULL ; v.split.temp <- NULL
              c.whole <- NULL ; c.split <- NULL 
              
              v.rep <- n-1
              
              #whole
              v.whole <- ( wfl-1)%*%t( wfl-1)
              v.whole <- c( wfl-1,v.whole[upper.tri(v.whole, diag=FALSE)])
              v.whole.denom <- prod( wfl)*n-1-sum(v.whole)-v.rep
              
              c.whole <- prod( wfl)*prod( sfl)*n/c( wfl, ( wfl%*%t( wfl))[upper.tri(( wfl)%*%t( wfl), diag=FALSE)])
              
              for (i in 1: length(v.whole)){
                pwr[n,i,deltao]<- (1-pf(qf((1-input$a4),v.whole[i], v.whole.denom),v.whole[i],v.whole.denom,ncp=ifelse(input$delta_type4==1,(Deltao[deltao]^2*(c.whole[i]*v.whole[i]))/(prod(sfl)+1),
                                                                                                                        (Deltao[deltao]^2*c.whole[i]/2)/(prod(sfl)+1))))#(Delta^2)*(c*nu1)

              }	
              
              #split
              v.split <- ( sfl-1)%*%t( sfl-1)
              v.split.temp <- (wfl-1)%*%t( sfl-1)
              v.split <- c( sfl-1, v.split[upper.tri(v.split, diag=FALSE)], as.vector(t(v.split.temp)))
              
              v.split.denom <- prod( wfl)*prod(sfl)*n-prod( wfl)*n-sum(v.split)
              
              c.split <- prod( wfl)*prod( sfl)*n/c( sfl, ( sfl%*%t( sfl))[upper.tri(( sfl)%*%t(sfl), diag=FALSE)], as.vector(t( wfl%*%t( sfl))))
              
              for (i in 1: length(v.split)){
                pwr[n,(i+length(v.whole)),deltao]<- (1-pf(qf((1-input$a4),v.split[i], v.split.denom),v.split[i],v.split.denom,ncp=ifelse(input$delta_type4==1,(Deltao[deltao]^2*(c.split[i]*v.split[i])) ,
                                                                                                                                          (Deltao[deltao]^2*c.split[i]/2) )))#(Delta^2)*(c*nu1)

              }
              
            }
            
            
          }
          
          
        }
        
        x<-c(1,2,3)
        j<-x[Deltao==input$plot_delta4]
        
        plot(2:100, pwr[2:100,1, j], type="l", ylim=c(0,1), xlim=c(0,min(100,n.choose+5)), ylab="Power", xlab="The number of replications", 
             main="The number of replications(r) vs Power",col=1, lwd=2)
        for (i in 2:ncol(Delta)) 
          lines(2:100, pwr[2:100,i,j], type="l", lty=i, lwd=2,col=i)
        
        abline(h=0.8, v=FF2$n,col="gray", lty=3)
        
        abline(h=0.9, col="gray", lty=3)
        legend("bottomright", legend=full_list, lty=seq(1:length(full_list)),col=seq(1:length(full_list)),lwd=2, adj=0)
        
      }
      
      
      )
##########################3
    }
    
    )
    
    
    
  }
