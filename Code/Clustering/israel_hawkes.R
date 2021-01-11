''' Data Preprocessing for Hawkes model of entire Mexico dataset'''

israel$date<-ymd(israel$EVENT_DATE)
israel$days<-yday(israel$date)
israel$totdays<-cumsum(israel$days)


startdate <- as.Date("2020-03-14", "%Y-%m-%d")

israel$num_days<-difftime(israel$date, startdate, units="days")
israel$num_days_duplicate<-duplicated(israel$num_days) | duplicated(israel$num_days, fromLast=TRUE)
israel$num_days_duplicate <- as.numeric(israel$num_days_duplicate == "TRUE")
israel$random <- ifelse(israel$num_days_duplicate == 0,0, israel$random <- runif(nrow(israel), min=0, max=1))
israel$day_rand <-israel$num_days+israel$random

israel<-israel[order(israel$day_rand),]
israel$day_rand<-as.numeric(israel$day_rand)

times_israel <- israel[,37]

# to export: write.csv(israel_events, "israel_events.csv", col.names=TRUE )

# hawkes model

pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_israel, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
israel_fit <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(israel_fit)

mfit_isr <- residuals(israel_fit, theoretical=TRUE, type="ord", std=TRUE, m= params(israel_fit)[1])

# residual/survival analysis
log.surv(mfit_isr)
#check stationarity: stationarity(mfit_isr, 10)


#calculate inter-arrival of residuals and check distribution
diff_res_isr <- diff(mfit_isr)
hist(diff_res_isr)
ks.exp.test(diff_res_isr, nrepl=1000)


''' data pre-processing for  cluster-based Hawkes models '''

mexico_c1 <- mexico_events[mexico_events$cluster==1,]
mexico_c2 <- mexico_events[mexico_events$cluster==2,]
mexico_c3 <- mexico_events[mexico_events$cluster==3,]
mexico_c4 <- mexico_events[mexico_events$cluster==4,]

startdate_m1 <- as.Date("2020-03-20", "%Y-%m-%d") # check that these start dates matches with the labels assigned to the clusters if you have generated yours
startdate_m2 <- as.Date("2020-03-20", "%Y-%m-%d")
startdate_m3 <- as.Date("2020-03-17", "%Y-%m-%d")
startdate_m4 <- as.Date("2020-03-17", "%Y-%m-%d")

# cluster 1

mexico_c1$date<-ymd(mexico_c1$EVENT_DATE)
mexico_c1$days<-yday(mexico_c1$date)
mexico_c1$totdays<-cumsum(mexico_c1$days)

mexico_c1$num_days<-difftime(mexico_c1$date, startdate_m1, units="days")
mexico_c1$num_days_duplicate<-duplicated(mexico_c1$num_days) | duplicated(mexico_c1$num_days, fromLast=TRUE)
mexico_c1$num_days_duplicate <- as.numeric(mexico_c1$num_days_duplicate == "TRUE")
mexico_c1$random <- ifelse(mexico_c1$num_days_duplicate == 0,0, mexico_c1$random <- runif(nrow(mexico_c1), min=0, max=1))
mexico_c1$day_rand <-mexico_c1$num_days+mexico_c1$random

mexico_c1<-mexico_c1[order(mexico_c1$day_rand),]
mexico_c1$day_rand<-as.numeric(mexico_c1$day_rand)

times_mexico_c1 <- mexico_c1[,12]


# cluster 2 
mexico_c2$date<-ymd(mexico_c2$EVENT_DATE)
mexico_c2$days<-yday(mexico_c2$date)
mexico_c2$totdays<-cumsum(mexico_c2$days)

mexico_c2$num_days<-difftime(mexico_c2$date, startdate_m2, units="days")
mexico_c2$num_days_duplicate<-duplicated(mexico_c2$num_days) | duplicated(mexico_c2$num_days, fromLast=TRUE)
mexico_c2$num_days_duplicate <- as.numeric(mexico_c2$num_days_duplicate == "TRUE")
mexico_c2$random <- ifelse(mexico_c2$num_days_duplicate == 0,0, mexico_c2$random <- runif(nrow(mexico_c2), min=0, max=1))
mexico_c2$day_rand <-mexico_c2$num_days+mexico_c2$random

mexico_c2<-mexico_c2[order(mexico_c2$day_rand),]
mexico_c2$day_rand<-as.numeric(mexico_c2$day_rand)

times_mexico_c2 <- mexico_c2[,12]

# cluster 3

mexico_c3$date<-ymd(mexico_c3$EVENT_DATE)
mexico_c3$days<-yday(mexico_c3$date)
mexico_c3$totdays<-cumsum(mexico_c3$days)

mexico_c3$num_days<-difftime(mexico_c3$date, startdate_m3, units="days")
mexico_c3$num_days_duplicate<-duplicated(mexico_c3$num_days) | duplicated(mexico_c3$num_days, fromLast=TRUE)
mexico_c3$num_days_duplicate <- as.numeric(mexico_c3$num_days_duplicate == "TRUE")
mexico_c3$random <- ifelse(mexico_c3$num_days_duplicate == 0,0, mexico_c3$random <- runif(nrow(mexico_c3), min=0, max=1))
mexico_c3$day_rand <-mexico_c3$num_days+mexico_c3$random

mexico_c3<-mexico_c3[order(mexico_c3$day_rand),]
mexico_c3$day_rand<-as.numeric(mexico_c3$day_rand)

times_mexico_c3 <- mexico_c3[,12]



# cluster 4

mexico_c4$date<-ymd(mexico_c4$EVENT_DATE)
mexico_c4$days<-yday(mexico_c4$date)
mexico_c4$totdays<-cumsum(mexico_c4$days)

mexico_c4$num_days<-difftime(mexico_c4$date, startdate_m4, units="days")
mexico_c4$num_days_duplicate<-duplicated(mexico_c4$num_days) | duplicated(mexico_c4$num_days, fromLast=TRUE)
mexico_c4$num_days_duplicate <- as.numeric(mexico_c4$num_days_duplicate == "TRUE")
mexico_c4$random <- ifelse(mexico_c4$num_days_duplicate == 0,0, mexico_c4$random <- runif(nrow(mexico_c4), min=0, max=1))
mexico_c4$day_rand <-mexico_c4$num_days+mexico_c4$random

mexico_c4<-mexico_c4[order(mexico_c4$day_rand),]
mexico_c4$day_rand<-as.numeric(mexico_c4$day_rand)

times_mexico_c4 <- mexico_c4[,12]



''' Cluster-based Hawkes models '''

# hawkes cluster 1

pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_mexico_c1, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
mexico_fit_c1 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(mexico_fit_c1)

mfit_mex1 <- residuals(mexico_fit_c1, theoretical=TRUE, type="ord", std=TRUE, m= params(mexico_fit_c1)[1])
log.surv(mfit_mex1)
stationarity(mfit_mex1, 10)


diff_res_mex1 <- diff(mfit_mex1)
hist(diff_res_mex1)
ks.exp.test(diff_res_mex1, nrepl=1000)

# hawkes cluster 2

pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_mexico_c2, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
mexico_fit_c2 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(mexico_fit_c2)

mfit_mex2 <- residuals(mexico_fit_c2, theoretical=TRUE, type="ord", std=TRUE, m= params(mexico_fit_c2)[1])
log.surv(mfit_mex2)
stationarity(mfit_mex2, 10)

diff_res_mex2 <- diff(mfit_mex2)
hist(diff_res_mex2)
ks.exp.test(diff_res_mex2, nrepl=1000)

# hawkes cluster 3

pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_mexico_c3, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
mexico_fit_c3 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(mexico_fit_c3)

mfti_mex3 <- residuals(mexico_fit_c3, theoretical=TRUE, type="ord", std=TRUE, m= params(mexico_fit_c3)[1])
log.surv(mfti_mex3)
stationarity(mfti_mex3, 10)

diff_res_mex3 <- diff(mfti_mex3)
hist(diff_res_mex3)
ks.exp.test(diff_res_mex3, nrepl=1000)

# hawkes cluster 4
pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_mexico_c4, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
mexico_fit_c4 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(mexico_fit_c4)

mfit_mex4 <- residuals(mexico_fit_c4, theoretical=TRUE, type="ord", std=TRUE, m= params(mexico_fit_c4)[1])
log.surv(mfit_mex4)
stationarity(mfit_mex4, 10)

diff_res_mex4 <- diff(mfit_mex4)
hist(diff_res_mex4)
ks.exp.test(diff_res_mex4, nrepl=1000)
