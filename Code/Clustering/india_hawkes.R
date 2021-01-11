''' Data Preprocessing for Hawkes model of entire Mexico dataset'''

##### ENTIRE DATASET

india$date<-ymd(india$EVENT_DATE)
india$days<-yday(india$date)
india$totdays<-cumsum(india$days)


startdate <- as.Date("2020-01-20", "%Y-%m-%d")

india$num_days<-difftime(india$date, startdate, units="days")
india$num_days_duplicate<-duplicated(india$num_days) | duplicated(india$num_days, fromLast=TRUE)
india$num_days_duplicate <- as.numeric(india$num_days_duplicate == "TRUE")
india$random <- ifelse(india$num_days_duplicate == 0,0, india$random <- runif(nrow(india), min=0, max=1))
india$day_rand <-india$num_days+india$random

india<-india[order(india$day_rand),]
india$day_rand<-as.numeric(india$day_rand)

times_india <- india[,37]

# to export: write.csv(india_events, "india_events.csv", col.names=TRUE )


# hawkes model
pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_india, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
india_fit <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(india_fit)

mfit_ind <- residuals(india_fit, theoretical=TRUE, type="ord", std=TRUE, m= params(india_fit)[1])

# residual/survival analysis
log.surv(mfit_ind)
# to check stationarity: stationarity(mfit_ind, 10)


# calculate inter-arrival of residuals and check distribution
diff_res_ind <- diff(mfit_ind)
hist(diff_res_ind)
ks.exp.test(diff_res_ind, nrepl=1000)



''' data pre-processing for  cluster-based Hawkes models '''

india_c1 <- india_events[india_events$cluster==1,]
india_c2 <- india_events[india_events$cluster==2,]
india_c3 <- india_events[india_events$cluster==3,]
india_c4 <- india_events[india_events$cluster==4,]

startdate_i1 <- as.Date("2020-02-01", "%Y-%m-%d") # check that these start dates match with the labels assigned to the clusters if you have generated yours
startdate_i2 <- as.Date("2020-01-20", "%Y-%m-%d")
startdate_i3 <- as.Date("2020-01-26", "%Y-%m-%d")
startdate_i4 <- as.Date("2020-02-14", "%Y-%m-%d")

# cluster 1

india_c1$date<-ymd(india_c1$EVENT_DATE)
india_c1$days<-yday(india_c1$date)
india_c1$totdays<-cumsum(india_c1$days)

india_c1$num_days<-difftime(india_c1$date, startdate_i1, units="days")
india_c1$num_days_duplicate<-duplicated(india_c1$num_days) | duplicated(india_c1$num_days, fromLast=TRUE)
india_c1$num_days_duplicate <- as.numeric(india_c1$num_days_duplicate == "TRUE")
india_c1$random <- ifelse(india_c1$num_days_duplicate == 0,0, india_c1$random <- runif(nrow(india_c1), min=0, max=1))
india_c1$day_rand <-india_c1$num_days+india_c1$random

india_c1<-india_c1[order(india_c1$day_rand),]
india_c1$day_rand<-as.numeric(india_c1$day_rand)

times_india_c1 <- india_c1[,12]


# cluster 2 
india_c2$date<-ymd(india_c2$EVENT_DATE)
india_c2$days<-yday(india_c2$date)
india_c2$totdays<-cumsum(india_c2$days)

india_c2$num_days<-difftime(india_c2$date, startdate_i2, units="days")
india_c2$num_days_duplicate<-duplicated(india_c2$num_days) | duplicated(india_c2$num_days, fromLast=TRUE)
india_c2$num_days_duplicate <- as.numeric(india_c2$num_days_duplicate == "TRUE")
india_c2$random <- ifelse(india_c2$num_days_duplicate == 0,0, india_c2$random <- runif(nrow(india_c2), min=0, max=1))
india_c2$day_rand <-india_c2$num_days+india_c2$random

india_c2<-india_c2[order(india_c2$day_rand),]
india_c2$day_rand<-as.numeric(india_c2$day_rand)

times_india_c2 <- india_c2[,12]

# cluster 3

india_c3$date<-ymd(india_c3$EVENT_DATE)
india_c3$days<-yday(india_c3$date)
india_c3$totdays<-cumsum(india_c3$days)

india_c3$num_days<-difftime(india_c3$date, startdate_i3, units="days")
india_c3$num_days_duplicate<-duplicated(india_c3$num_days) | duplicated(india_c3$num_days, fromLast=TRUE)
india_c3$num_days_duplicate <- as.numeric(india_c3$num_days_duplicate == "TRUE")
india_c3$random <- ifelse(india_c3$num_days_duplicate == 0,0, india_c3$random <- runif(nrow(india_c3), min=0, max=1))
india_c3$day_rand <-india_c3$num_days+india_c3$random

india_c3<-india_c3[order(india_c3$day_rand),]
india_c3$day_rand<-as.numeric(india_c3$day_rand)

times_india_c3 <- india_c3[,12]



# cluster 4

india_c4$date<-ymd(india_c4$EVENT_DATE)
india_c4$days<-yday(india_c4$date)
india_c4$totdays<-cumsum(india_c4$days)

india_c4$num_days<-difftime(india_c4$date, startdate_i4, units="days")
india_c4$num_days_duplicate<-duplicated(india_c4$num_days) | duplicated(india_c4$num_days, fromLast=TRUE)
india_c4$num_days_duplicate <- as.numeric(india_c4$num_days_duplicate == "TRUE")
india_c4$random <- ifelse(india_c4$num_days_duplicate == 0,0, india_c4$random <- runif(nrow(india_c4), min=0, max=1))
india_c4$day_rand <-india_c4$num_days+india_c4$random

india_c4<-india_c4[order(india_c4$day_rand),]
india_c4$day_rand<-as.numeric(india_c4$day_rand)

times_india_c4 <- india_c4[,12]





''' Cluster-based Hawkes models '''

# hawkes cluster 1

pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_india_c1, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
india_fit_c1 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(india_fit_c1)

mfit_ind1 <- residuals(india_fit_c1, theoretical=TRUE, type="ord", std=TRUE, m= params(india_fit_c1)[1])
log.surv(mfit_ind1)
stationarity(mfit_ind1, 10)

diff_res_ind1 <- diff(mfit_ind1)
hist(diff_res_ind1)
ks.exp.test(diff_res_ind1, nrepl=1000)

# hawkes cluster 2

pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_india_c2, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
india_fit_c2 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(india_fit_c2)

mfit_ind2 <- residuals(india_fit_c2, theoretical=TRUE, type="ord", std=TRUE, m= params(india_fit_c2)[1])
log.surv(mfit_ind2)
stationarity(mfit_ind2, 10)

diff_res_ind2 <- diff(mfit_ind2)
hist(diff_res_ind2)
ks.exp.test(diff_res_ind2, nrepl=1000)

# hawkes cluster 3

pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_india_c3, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
india_fit_c3 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(india_fit_c3)

mfti_ind3 <- residuals(india_fit_c3, theoretical=TRUE, type="ord", std=TRUE, m= params(india_fit_c3)[1])
log.surv(mfti_ind3)
stationarity(mfti_ind3, 10)

diff_res_ind3 <- diff(mfti_ind3)
hist(diff_res_ind3)
ks.exp.test(diff_res_ind3, nrepl=1000)


# hawkes cluster 4
pstart = c(mu = .004, C = .004, a = 0.5)
ppm <- ptproc(pts = times_india_c4, cond.int = hawkes.cond.int, params = pstart)
condition(ppm) <- penalty(code = NULL, condition = "any(params < 0)")
india_fit_c4 <- ptproc.fit(ppm, optim.control = list(trace = 2), alpha = 1e+5, hessian = TRUE)
summary(india_fit_c4)

mfit_ind4 <- residuals(india_fit_c4, theoretical=TRUE, type="ord", std=TRUE, m= params(india_fit_c4)[1])
log.surv(mfit_ind4)
stationarity(mfit_ind4, 10)

diff_res4_ind <- diff(mfit_ind4)
hist(diff_res4_ind)
ks.exp.test(diff_res4_ind, nrepl=1000)
