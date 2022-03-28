# Plot moonrise & moonset relative to sunrise & sunset each day at each site
library(suncalc)


start = '2016-05-01 00:00:00'
end = '2019-04-30 23:59:59'
lats = as.numeric(c(41.06165,40.22999,39.83295,
                    39.19192,38.37337,37.16452,
                    35.30183,33.66992,32.10527,
                    30.58295,30.27818))
lons = as.numeric(c(-66.35155,-67.97798,-69.98194,
                    -72.22735,-73.36985,-74.46585,
                    -74.87895,-75.9977,-77.09067,
                    -77.39002,-80.22085))


lunTime = getMoonTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd+(60*60*24)),by=1),
                       lat=lats[j],lon=lons[j],
                       keep=c("rise","set"),
                       tz="UTC")

dayData = getSunlightTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd),by=1),
                           lat=lats[j],lon=lons[j],
                           keep=c("sunrise","sunset"),
                           tz="UTC")

