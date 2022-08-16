#Name: Cox_model_function.R----
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Fit an arbitrarily large number of predictors (genes) to a
#cox model to predict survival time

#cox_model_fitter----
cox_model_fitter <- function(my.seed       =1,
                             cox.df        =NULL,
                             gene.num      =1800,
                             cox.predictors=NULL,
                             tumor.stage   =FALSE,
                             tumor.n       =FALSE,
                             tumor.m       =FALSE,
                             regular.cox   =TRUE,
                             save.regular.cox.genes =TRUE,
                             remove.stage   = c("tumor.stage1","tumor.stage2","tumor.stage3", "tumor.stage4"),
                             remove.n.stage = c("ajcc.n0", "ajcc.n1", "ajcc.n2", "ajcc.n3"),
                             my.filename   ="my_saved_genes.csv"){
  
  #Doing input sanity checks----
  if(missing(my.seed)){
    stop("You must specify a seed for reproducible analysis.")
  }
  
  if(missing(cox.df)){
    stop("You must include a data.frame to carry out cox analysis.")
  }
  
  if(missing(cox.predictors)){
    stop("You must include predictors for your cox model.")
  }
  
  if (class(my.seed)!="numeric"){
    stop("You must specify a number for your seed.")
  }
  
  if(class(cox.df)!="data.frame"){
    stop("Your cox.df parameter must be of type 'data.frame'.")
  }
  
  #If the packages are installed, they----
  #will be loaded. If they are not,
  #the packages will be installed
  #from CRAN and then loaded.
  require(BiocGenerics)
  require(glmnet)
  require(lmtest)
  require(survival)
  require(survminer)
  
  #Setting the seed for reproducible----
  #output
  set.seed(my.seed)
  
  #Making the list that will store the----
  #data we will return
  cox_data <- list()
  
  #The predictors for the cox model----
  my_predictors <- cox.predictors
  
  my_predictors <- head(my_predictors, n=gene.num)
  my_predictors <- sapply(my_predictors, gsub, pattern="-",replacement=".")
  my_predictors <- unlist(my_predictors)
  colname_changes <- sapply(colnames(cox.df), gsub, pattern="-",replacement=".")
  colname_changes <- sapply(colnames(cox.df), gsub, pattern="_",replacement=".")
  colname_changes <- sapply(colnames(cox.df), gsub, pattern="/",replacement=".")
  colname_changes <- unlist(colname_changes)
  colnames(cox.df) <- colname_changes
  my_predictors <- intersect(my_predictors, colnames(cox.df))
  my_predictors <- paste("~", paste(my_predictors[1:length(my_predictors)], collapse = "+"))
  
  
  if(tumor.stage==TRUE & tumor.n==FALSE & tumor.m==FALSE){
    my_predictors <- paste(my_predictors, "tumor.stage", sep = "+")
    my_predictors <- as.formula(my_predictors)
  }else if (tumor.stage==TRUE & tumor.n==TRUE & tumor.m==FALSE){
    my_predictors <- paste(my_predictors,"ajcc.n", sep = "+")
    my_predictors <- as.formula(my_predictors)
  }else if (tumor.stage==TRUE & tumor.n==TRUE & tumor.m==TRUE){
    my_predictors <- paste(my_predictors, "ajcc.m", sep = "+")
    my_predictors <- as.formula(my_predictors)
  }else{
    my_predictors <- as.formula(my_predictors)
    # print("This is the genes only predictor....")
  }
  
  my_x <- model.matrix(my_predictors, cox.df)
  #The response object for the cox model----
  my_y <- Surv(time = cox.df$days.to.last.follow.up, event = cox.df$vital.status)
  #The 10-fold cross-validation fit----
  cv_fit <- cv.glmnet(x = my_x, y = my_y, nfolds = 10, type.measure = "C", maxit=100000, family="cox", parallel = TRUE)
  
  
  #Looking to see which genes are the most important
  fit <- glmnet(my_x, my_y, family =  "cox", maxit = 100000)
  Coefficients <- coef(fit, s = cv_fit$lambda.min)
  
  Active.Index <- which(as.logical(Coefficients) != 0)
  Active.Coefficients  <- Coefficients[Active.Index]
  active_genes <-rownames(Coefficients)[Active.Index]
  
  if(tumor.stage==TRUE & tumor.n==FALSE & tumor.m==FALSE){
    print("This is the genes + tumor stage predictor...")
    active_genes <- active_genes[!active_genes %in% remove.stage]
    active_genes <- c(active_genes, "tumor.stage")
  }
  
  if(tumor.stage==TRUE & tumor.n==TRUE & tumor.m==FALSE){
    print("This is the genes + tumor stage + n stage predictor...")
    active_genes <- active_genes[!active_genes %in% remove.n.stage]
    active_genes <- c(active_genes, "tumor.stage", "ajcc.n")
  }
  
  #Getting survival
  my_surv <- survival::survfit(cv_fit, s = "lambda.min", x = my_x, y = my_y)
  
  #Adding the relevant data bits to list to return
  cox_data[["CV"]] <- cv_fit
  cox_data[["Coefficients"]] <- Coefficients
  cox_data[["Active Coefficients"]] <- Active.Coefficients
  cox_data[["Active Index"]] <- Active.Index
  cox_data[["Active Genes"]] <- active_genes
  cox_data[["Predictors"]] <- my_predictors
  cox_data[["Predicted Survival"]] <- my_surv
  
  
  #Regular Cox
  if(regular.cox==TRUE){
    active_predictors <-paste(active_genes, collapse = "+")
    regular_cox_df <- cox.df[,active_genes]
    regular_cox_df$days.to.last.follow.up <- cox.df$days.to.last.follow.up
    regular_cox_df$vital.status <- cox.df$vital.status
    my_formula <- paste("~", paste(active_genes[1:length(active_genes)], collapse = "+"))
    regular_cox <- coxph(Surv(time = regular_cox_df$days.to.last.follow.up, event = regular_cox_df$vital.status)~., data = regular_cox_df)
    cox_data[["Regular Cox"]] <- regular_cox
    cox_data[["Summary of Regular Cox"]] <- summary(regular_cox)
    if(save.regular.cox.genes==TRUE){
      summ_regular_cox <- summary(regular_cox)
      write.csv(summ_regular_cox$coefficients, file = my.filename)
    }
    
  }
  
  #Returning our finished output----
  return(cox_data)

}



