### Data Preparation function
dataprep = function(filename, age, relative_path = 'data/'){
    #######################################
    ### DATA for PCA
    #######################################
    ## Read in raw data
    dat.all = read.csv(paste0(relative_path, filename, ".csv"), fileEncoding = 'UTF-8-BOM')
    dat.all = as.data.frame(dat.all)
    ## Exclude GMCSF and Exclude the covariates with insufficient observations
    ind = match(c("GMCSF", "CCL4", "CCL17"),colnames(dat.all)); ind = ind[which(!is.na(ind))]; ind
    if(length(ind)>0) {dat.all = dat.all[, -ind]}; head(dat.all)
    
    ## Retrieve Response names
    temp = read.csv(paste0(relative_path, filename, ".csv"), fileEncoding = 'UTF-8-BOM',check.names=FALSE) 
    ind = match(c("GMCSF", "CCL4", "CCL17"),colnames(temp)); ind = ind[which(!is.na(ind))]; ind
    if(length(ind)>0) {temp = temp[, -ind]}; head(temp)
    
    ## Retrieve header
    ind.header = which(colnames(dat.all)%in%c("Sex", "Mouse.ID","Treatment.ID" ,"Treatment", "Age", "Date"))
    map = read.csv(paste0(relative_path,"Response_map.csv"), fileEncoding = 'UTF-8-BOM'); 
    resp.map = data.frame(
      Response = colnames(dat.all)[-ind.header],
      ResponsePlot = gsub("cDC", "DC",gsub("[|].*","",gsub("\\s", "",colnames(temp))))[-ind.header]) %>% 
      left_join(map, by = "ResponsePlot") %>%
      dplyr::select("Response.ID", "Response","ResponsePlot");resp.map
    resp.map = resp.map[order(resp.map$Response.ID),]; resp.map
    resp.ls = colnames(dat.all)[-ind.header]; resp.ls
    
    # Retrieve Treatment names
    map = read.csv(paste0(relative_path,"Treatment_map.csv"), fileEncoding = 'UTF-8-BOM')
    trt.map = map %>% dplyr::filter(Treatment.ID %in% dat.all$Treatment.ID); trt.map
    dat.all$Treatment.ID = factor(dat.all$Treatment.ID)
    dat.all$Mouse.ID = factor(dat.all$Mouse.ID)
    if(0 %in% dat.all$Age){
      dat.all$Age[which(dat.all$Age == 0)] =  "Young"
    }
    if(1 %in% dat.all$Age){
      dat.all$Age[which(dat.all$Age == 1)] =  "Aged"
    }
    dat.all$Age = factor(dat.all$Age, levels = c("Aged", "Young"))
    
    # Check point for raw dataset
    unique(dat.all$Mouse.ID)
    unique(dat.all$Age)
    dat.all %>% dplyr::select(Mouse.ID, Age)
    table(dat.all$Mouse.ID)
    table(dat.all$Age)
    table(dat.all$Treatment.ID)
    identical(sapply( dat.all$Mouse.ID, function(x){grepl("A", x, fixed = TRUE)}),(dat.all$Age == "Aged"))
    
    ind.header = which(colnames(dat.all)%in%c("Sex", "Mouse.ID","Treatment.ID" ,"Treatment", "Age", "Date"))
    res.fit = list()
    
    #######################################
    ### DATA for SINGLE-AGE ANALYSIS
    #######################################
    if(age %in% c("Aged", "Young")){
    dat1 = dat.all %>% dplyr::filter(Age == age)
    dat1[,-ind.header] = apply(dat1[,-ind.header], 2, as.numeric)
    # #### Data for plotting  increment/decrement for each group
    # res.all2 = c()
    # for(ycov in resp.ls){
    #   # Model fitting
    #   print(paste0("Now testing for ",ycov))
    #   o = lmer(paste0("log(",ycov,"+1)", "~ (1|Mouse.ID) +
    #                             Treatment.ID"),data = dat1)
    #   anova(o, ddf = "Kenward-Roger") # same result in SAS
    #   res2  = c()
    #   for (ind.trt in trt.map[-1,1]){
    #     temp = rep(0,nrow(trt.map)-1)
    #     temp[which(ind.trt == trt.map)-1] = 1
    #     res2 = rbind(res2,
    #                  # Aged(Trt)
    #                  contest(model = o,
    #                          L = c(1, # intercept
    #                                # 0, # Young -Aged
    #                                temp # Treatment2.. - Unsti
    #                                # rep(0, nrow(trt.map)-1) # Young:Treatment2..-Aged:Treatment2..
    #                          ),
    #                          joint = FALSE, # F for single t test, T for joint F test
    #                          confint = TRUE,
    #                          ddf = "Kenward-Roger"),
    #                  #Aged(Unsti)
    #                  contest(model = o,
    #                          L = c(1, # intercept
    #                                # 0, # Young -Aged
    #                                rep(0, nrow(trt.map)-1) # Treatment2.. - Unsti
    #                                # rep(0, nrow(trt.map)-1) # Young:Treatment2..-Aged:Treatment2..
    #                          ),
    #                          joint = FALSE, # F for single t test, T for joint F test
    #                          confint = TRUE,
    #                          ddf = "Kenward-Roger"))
    #   }
    #   res = data.frame(Estimate = res2$Estimate)
    #   res$Items =factor(rep(c("Trt", "Unsti"), times = nrow(trt.map)-1), levels = c("Unsti","Trt"))
    #   res$Age = rep(c(age,age), times = nrow(trt.map)-1)
    #   res$Response = resp.map[which(ycov == resp.map[,2]),3]
    #   res$Response.ID = resp.map[which(ycov == resp.map[,2]),1]
    #   res$Treatment.ID = rep(trt.map[-1,1], each = 2)
    #   res$Treatment = rep(trt.map[-1,2],each = 2)
    #   res.all2 = rbind(res.all2, res)
    # }
    # save(res.all2, file = paste0("test.res/", filename, age,".res.all2.rda"))
    
    #### Data for plotting test of Treatment effect
    res.all3 = c()
    for(ycov in resp.ls){
      # Model fitting
      print(paste0("Now testing for ",ycov))
      # pure linear
      # o = lm(paste0("log(",ycov,"+1)", "~
      #                           Treatment.ID "),data = dat1)

      # mixed linear
      # o = lmer(paste0("log(",ycov,"+1)", "~ (1|Mouse.ID) + (1|Treatment.ID) +
      #                           Treatment.ID "),data = dat1)
      o = lmer(paste0("log(",ycov,"+1)", "~ (1|Mouse.ID) +
                                Treatment.ID "),data = dat1)
      ind.obs = which(!is.na(dat1%>%dplyr::select(all_of(ycov))))
      # anova(o, ddf = "Kenward-Roger") # same result in SAS
      # Check residuals
      plot(predict(o,re.form = NULL) ,predict(o,re.form = NA) , main = ycov)
      tempplot = data.frame(fitted = predict(o, re.form = NA), residual = resid(o, type = "response"), Mouse.ID= dat1$Mouse.ID[ind.obs], Treatment.ID= dat1$Treatment.ID[ind.obs])
      fig = ggplot(tempplot, aes(fitted,residual, shape = Mouse.ID, colour = Treatment.ID))+
        geom_point(aes(group = interaction(Mouse.ID, Treatment.ID)))
      print(fig)
      # par(oma=c(0, 0, 0, 5))
      # plot(x = fitted(o), y = resid(o, type = "response"), main = ycov, col = factor(dat1$Mouse.ID), pch = factor(dat1$Treatment.ID), ylab = "conditional residuals")
      # legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
      #        legend = levels(factor(dat1$Mouse.ID)), col = 1:20, cex = 0.8)
      # legend("bottomright", legend = levels(factor(dat1$Treatment.ID)), pch = 1:20, cex = 0.8)
      
      # next
      
    
      res.fit[[ycov]] = o
      res3  = c()
      for (ind.trt in trt.map[-1,1]){
        temp = rep(0,nrow(trt.map)-1)
        temp[which(ind.trt == trt.map)-1] = 1
        res3 = rbind(res3,
                     # Aged(Trt) - Aged(Unsti)
                     contest(model = o,
                             L = c(0, # intercept
                                   # 0, # Young -Aged
                                   temp# Treatment2.. - Unsti
                                   # rep(0, nrow(trt.map)-1) # Young:Treatment2..-Aged:Treatment2..
                             ),
                             joint = FALSE, # F for single t test, T for joint F test
                             confint = TRUE,
                             ddf = "Kenward-Roger"))
      }
      res = (res3)
      res$Items =factor(rep(paste0("Trt-Unsti|", age), times = nrow(trt.map)-1))
      res$Response = resp.map[which(ycov == resp.map[,2]),3]
      res$Response.ID = resp.map[which(ycov == resp.map[,2]),1]
      res$Treatment.ID = rep(trt.map[-1,1], each = 1)
      res$Treatment = rep(trt.map[-1,2],each = 1)
      res.all3 = rbind(res.all3, res)
    }
    save(res.all3, file = paste0("test.res/", filename,age,".res.all3.rda"))
    save(res.fit, file = paste0("test.res/", filename,age,".res.fit.M1.rda"))
    }
      
    #######################################
    ### DATA for BOTH-AGE ANALYSIS
    #######################################
    if(age == "AgedYoung"){
      res.all = c()
      dat.all[,-ind.header] = apply(dat.all[,-ind.header], 2, as.numeric)
      for(ycov in resp.ls){# "MIP1b", "TARC" Exclude the covariates with insufficient observations
      #### Data for plotting test of Treatment effect
      # Model fitting
      print(paste0("Now testing for ",ycov))
      # o = lm(paste0("log(", ycov, "+1) ~ Age  + Treatment.ID + Treatment.ID:Age"),data = dat.all)
      o = lmer(paste0("log(", ycov, "+1) ~ Age  + Treatment.ID + Treatment.ID:Age +  (1|Mouse.ID) "),data = dat.all)
      
      res.fit[[ycov]] = o
      ind.obs = which(!is.na(dat.all%>%dplyr::select(all_of(ycov))))
      
      # anova(o, ddf = "Kenward-Roger") # same result in SAS
      # Check residuals
      plot(predict(o,re.form = NULL) ,predict(o,re.form = NA) , main = ycov)
      tempplot = data.frame(fitted = predict(o, re.form = NA), residual = resid(o, type = "response"), Mouse.ID= dat.all$Mouse.ID[ind.obs], Treatment.ID= dat.all$Treatment.ID[ind.obs])
      fig = ggplot(tempplot, aes(fitted,residual, shape = Mouse.ID, colour = Treatment.ID))+
        geom_point(aes(group = interaction(Mouse.ID, Treatment.ID)))
      print(fig)
      
      # par(oma=c(0, 0, 0, 5))
      # plot(x = fitted(o), y = resid(o, type = "response"), main = ycov, col = factor(dat.all$Age), ylab = "conditional residuals")
      # legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
      #        legend = levels(factor(dat.all$Age)), col = 1:20, cex = 0.8, pch = 1:20)
      # 
      # next
      summary(o)
      ls_means(o)
      ranef(o)
      fixef(o)
      print(summary(o)$sigma)
      # Testing Constrasts: Young(Unsti)-Aged(Unsti)
      res1 = contest(model = o,
                     L = c(0, # intercept
                           1, # Young -Aged
                           rep(0, nrow(trt.map)-1),# Treatment2.. - Unsti
                           rep(0, nrow(trt.map)-1) # Young:Treatment2..-Aged:Treatment2..
                     ),
                     joint = FALSE, # F for single t test, T for joint F test
                     confint = TRUE,
                     ddf = "Kenward-Roger")
      
      # Testing Constrasts: Young(Trt)-Aged(Trt)
      res2  = c()
      for (ind.trt in trt.map[-1,1]){
        temp = rep(0,nrow(trt.map)-1)
        temp[which(ind.trt == trt.map)-1] = 1
        # Young(Trt1)-Aged(Trt1)
        res2 = rbind(res2,
                     contest(model = o,
                             L = c(0, # intercept
                                   1, # Young -Aged
                                   rep(0, nrow(trt.map)-1), # Trt.. - Unsti
                                   temp  # Young:Trt-Aged:Trt
                             ),
                             joint = FALSE, # F for single t test, T for joint F test
                             confint = TRUE,
                             ddf = "Kenward-Roger"),
                     # Testing Constrasts:Young(Trt-Unsti) - Aged(Trt -Unsti)
                     contest(model = o,
                             L = c(0, # intercept
                                   0, # Young -Aged
                                   rep(0, nrow(trt.map)-1), # Treatment2.. - Unsti
                                   temp # Young:Treatment2..-Aged:Treatment2..
                             ),
                             joint = FALSE, # F for single t test, T for joint F test
                             confint = TRUE,
                             ddf = "Kenward-Roger"),
                     # Testing Constrasts: Aged(Trt) - Aged(Unsti)
                     contest(model = o,
                             L = c(0, # intercept
                                   0, # Young -Aged
                                   temp, # Treatment2.. - Unsti
                                   rep(0, nrow(trt.map)-1) # Young:Treatment2..-Aged:Treatment2..
                             ),
                             joint = FALSE, 
                             confint = TRUE,
                             ddf = "Kenward-Roger"),
                     # Testing Constrasts: Young(Trt) - Young(Unsti)
                     contest(model = o,
                             L = c(0, # intercept: Aged
                                   0, # Young -Aged
                                   temp, # Treatment2.. - Unsti
                                   temp  # Young:Treatment2..-Aged:Treatment2..
                             ),
                             joint = FALSE, # F for single t test, T for joint F test
                             confint = TRUE,
                             ddf = "Kenward-Roger"))
      }
      res = rbind(res1, res2)
      res$Response =  resp.map[which(ycov == resp.map[,2]),3]
      res$Response.ID = resp.map[which(ycov == resp.map[,2]),1]
      res$Treatment.ID = c("Unsti", rep(trt.map[-1,1],each = 4))
      res$Treatment = c("Unsti", rep(trt.map[-1,2],each = 4))
      res$Items = c(paste0("AgeDiff"), rep(c("AgeDiff|Trt", "AgeReactDiff|Trt","Trt-Unsti|Aged","Trt-Unsti|Young"),times =nrow(trt.map)-1))
      res$YoungHigher = res$Estimate>0
      res$Significance = res$`Pr(>|t|)` <0.05
      # print(res)
      res.all = rbind(res.all, res)
      }
      save(res.all, file = paste0("test.res/", filename,age,".res.all.rda"))
      save(res.fit, file = paste0("test.res/", filename,age,".res.fit.M2.rda"))
    }
    

}

