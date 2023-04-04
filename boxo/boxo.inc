       integer*4 ii,jj                                                          
       real*8 time                                                              
       real*8 k1mean                                                            
       real*8 e1                                                                
       real*8 k2mean                                                            
       real*8 e2                                                                
       real*8 tmean                                                             
       integer*4 case                                                           
       real*8 temp                                                              
c                                                                               
c                                                                               
       time   = xdata(1,iobs)                                                   
       k1mean = gpar(1)                                                         
       e1     = gpar(2)                                                         
       k2mean = gpar(3)                                                         
       e2     = gpar(4)                                                         
       tmean  = gpar(5)                                                         
       case   = gpar(6)                                                         
       temp   = lpar(1)                                                         
