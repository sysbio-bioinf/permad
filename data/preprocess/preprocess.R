## preprocessing of PERMAD data cohort 1 and cohort 2

## load raw data from Myriad
library(openxlsx)
raw.1 <- read.xlsx("PERMAD_cohort_1.xlsx", sheet=1, startRow=8, detectDates=T)
raw.2 <- read.xlsx("PERMAD_cohort_2.xlsx", sheet=1, startRow=8)

## extract feature names
featurenames <- names(read.xlsx("PERMAD_cohort_1.xlsx", sheet=1, colNames=T, rows=1, sep.names=" "))[-1]

## some patients have to be removed, i.e. grey boxes in Myriad xlsx with "kein Progress ...", "Einwilligung ...", "nicht auswertbar"
dropouts.1.ID <- c("106-001", "110-004", "112-001", "204-001")
dropouts.2.ID <- c("107-010", "109-004", "112-002", "113-002", "113-005", "201-010")
dropouts.1.myriadID <- sapply(strsplit(raw.1[match(dropouts.1.ID, raw.1$ID), "Myriad.ID.+.Samples"], " "), "[",1)
dropouts.2.myriadID <- sapply(strsplit(raw.2[match(dropouts.2.ID, raw.2$ID), "Myriad.ID.+.Samples"], " "), "[",1)
dropouts.1 <- which(!is.na(sapply(raw.1[, "Myriad.ID.+.Samples"], \(x) match(strsplit(x, " ")[[1]][1], dropouts.1.myriadID))))
dropouts.2 <- which(!is.na(sapply(raw.2[, "Myriad.ID.+.Samples"], \(x) match(strsplit(x, " ")[[1]][1], dropouts.2.myriadID))))

## last rows have to be removed, i.e. "Comments ..."
droprows.1 <- which(is.na(raw.1[,1]))
droprows.2 <- which(is.na(raw.2[,1]))

## create data matrix, i.e. combine data from cohort 1 and 2, without dropouts and droprows
dat.1 <- t(raw.1[-c(dropouts.1, droprows.1), -c(1:5)])
dat.2 <- t(raw.2[-c(dropouts.2,droprows.2), -c(1:5)])

dat <- cbind(dat.1, dat.2)
rownames(dat) <- featurenames
colnames(dat) <- sub(" ", "", sub(" ", "_", c(raw.1[-c(dropouts.1, droprows.1), "Myriad.ID.+.Samples"], raw.2[-c(dropouts.2, droprows.2), "Myriad.ID.+.Samples"])))

## clean data: 1) set signals below/above detection limit to LLOQ/max value from Myriad 2) remove spaces 3) convert to numeric
dat.cleaned <- sub("<", "", dat)
dat.cleaned <- sub(">", "", dat.cleaned)
dat.cleaned <- gsub(" ", "", dat.cleaned)
dat.cleaned <- apply(dat.cleaned, 1:2, as.numeric)

## collect info about cohort and sample id
cohort <- c(rep(1, ncol(dat.1)), rep(2, ncol(dat.2)))
sample.id <- as.numeric(c(sapply(strsplit(raw.1[-c(dropouts.1, droprows.1), "Myriad.ID.+.Samples"], " "), "[", 1), sapply(strsplit(raw.2[-c(dropouts.2, droprows.2), "Myriad.ID.+.Samples"], " "), "[", 1)))

## extract info about before treatment, ct, last ct (i.e. green/yellow/red in Myriad xlsx)
time.deltact <- as.numeric(c(raw.1$Abstand.CT.in.Tagen[-c(dropouts.1,droprows.1)], raw.2$Abstand.CT.in.Tagen[-c(dropouts.2,droprows.2)]))
time.deltact[max(which(sample.id==26))] <- 15 # time for progress for 26 21 was 22.06.2016, i.e. +15, see "Ettrich permad sample uebersicht.txt".
trtct <- (!is.na(time.deltact))*2 # ct
trtct[!duplicated(sample.id)] <- 1 # before treatment
btonlyone <- c(9, 10, 12, 15, 32, 33, 35, 36, 38, 39, 40, 46, 49, 50)
trtct[which(!duplicated(sample.id)&!(sample.id %in% btonlyone))+1] <- 1 # second also before treatment
trtct[c(which(!duplicated(sample.id))-1,length(trtct))[-1]] <- 3 # mark last ct as special

## normalize, i.e. substract mean of before treatment samples (green in Myriad xlsx)
dat.normalized <- dat.cleaned
for(i in unique(sample.id))
    dat.normalized[, sample.id==i] <- dat.cleaned[, sample.id==i] - rowMeans(dat.cleaned[, sample.id==i & trtct==1, drop=F])

## create different time vectors from date, e.g. time before progress (bp), time of ct before progress (ct), timestamp, ...
time.1 <- as.Date(raw.1$Datum.Abnahme[-c(dropouts.1, droprows.1)], format="%d/%m/%Y")
time.1[is.na(time.1)] <- raw.1$Datum.Abnahme[-c(dropouts.1, droprows.1)][is.na(time.1)]
time.2 <- as.Date(raw.2$Datum.Abnahme[-c(dropouts.2, droprows.2)], format="%d/%m/%Y")
time.date <- c(time.1, time.2)

time.bp <- unsplit(mapply(\(x, y, z) difftime(x, x[y==3], units="days") - z[y==3], split(time.date, sample.id), split(trtct, sample.id), split(time.deltact, sample.id)), sample.id)
time.ct <- unsplit(mapply(\(x, y, z) difftime(x, x[y==3], units="days") - z[y==3] + z, split(time.date, sample.id), split(trtct, sample.id), split(time.deltact, sample.id)), sample.id)
time.id <- unsplit(lapply(split(time.bp, sample.id), \(x) rev(abs(seq.int(-46, along.with=x)))), sample.id)

## labels, i.e. split at 100
labs <- (time.bp > -100) +1

## we have progress ct samples of 3 patients (4 13, 10 8, 21 7) with time of sample earlier than ct. After discussion with Ettrich we have to remove two of them, see "Ettrich permad sample uebersicht.txt".
## we also remove all before treatment samples, as we don't need them any more
invalid <- match(c("10_8", "21_7"), colnames(dat.normalized))
remove <- c(invalid, which(trtct==1))

## create large RDATA object with data, labs, sample.id, cohort, time.id, time, time.ct, features, name
PERMAD <- list(data=dat.normalized[,-remove], labs=labs[-remove], sample.id=sample.id[-remove], cohort=cohort[-remove], time=as.numeric(time.bp[-remove]), time.ct=as.numeric(time.ct[-remove]), time.id=as.numeric(time.id[-remove]), features=featurenames, name="PERMAD")
save(PERMAD, file="PERMAD.RData")
