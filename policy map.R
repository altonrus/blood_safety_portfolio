library(ggplot2)
library(data.table)
# library(mapproj)
library(scales)
library(sf)
library(rmapshaper)
library(maptools)
library(grid)
library(gtable)
# library(rgdal)
# library(tidyverse)
# library(geojsonio)
library(RColorBrewer)
theme_set(theme_minimal())
# library(cowplot)
# library(gridExtra)
# library(rgeos)
# library(dplyr)
# library(concaveman)
# usa_zip5 <- data.table(st_read('data/ZCTA_2010Census_DP1.shp'))
# dt <- usa_zip5[ , c("ZCTA5CE10", "geometry")]
# dt[, zip3 := substr(ZCTA5CE10, 1, 3)]
# 
# dt_zip3 <- dt[, .(geom = st_union(geometry)),by=zip3]
# test <- as(as_Spatial(dt_zip3$geom, IDs = dt_zip3$zip3), "SpatialPolygonsDataFrame" )
# writePolyShape(test, fn = "data/dt_zip3.shp")

zip3_shapes <- readShapePoly(fn = "data/dt_zip3.shp")
# 
# zip3_shapes_PR <- zip3_shapes[zip3_shapes$SP_ID %in% c("006", "007", "008"),]
# zip3_shapes_PR <- elide(zip3_shapes_PR, shift = c(-5, 5))
# 
# zip3_shapes_AK <- zip3_shapes[zip3_shapes$SP_ID %in% c("995", "996", "997", "998", "999"),]
# zip3_shapes_AK <- elide(zip3_shapes_PR, shift = c(0, -25))
# 
# zip3_shapes_HI <- zip3_shapes[zip3_shapes$SP_ID %in% c("967", "968"),]
# zip3_shapes_HI <- elide(zip3_shapes_PR, shift = c(5, 25))
# 
# zip3_shapes <- zip3_shapes[!(zip3_shapes$SP_ID %in% c("006", "007", "008",
#                                                       "967", "968",
#                                                       "995", "996", "997", "998", "999")),]
# 
# zip3_shapes <- rbind.SpatialPolygonsDataFrame(zip3_shapes,zip3_shapes_AK, zip3_shapes_HI, zip3_shapes_PR)


zip3_sf <- st_as_sf(zip3_shapes)

# Sys.time()
# zip3_sf <- ms_simplify(zip3_sf)
# Sys.time()

dt_zip3 <- data.table(
  zip3_sf
)

setnames(dt_zip3, "SP_ID", "zip3")
dt_zip3$dummy <- NULL



# #plot(dt_zip3)
# zip3_shapeS_simp <- ms_simplify(zip3_shapes, keep_shapes = TRUE)

# zip3_shapes <- left_join(zip3_shapes, psa_percents_zip[year == "2017" & season == "H", c("year", "a1", "a2", "a3", "a4", "zip3", "STUSPS")])

#dt <- data.table(data.frame(zip3 = zip3_shapes$SP_ID))



# Read in PSA percents
psa_percents_zip <- fread("results/PSA_zip3_opt_by_group_combined.csv")[ , lapply(.SD, sum), .SD = 3:8, by = c("group", "year")]
psa_percents_zip[ , z := (10000 - z)/10000]
psa_percents_zip[, c("a1", "a2", "a3", "a4", "m1") :=lapply(.SD, function(x) x/10000), .SDcols=4:8]
psa_percents_zip[ , STUSPS := substr(group, 1, 2)]
psa_percents_zip[ , season := substr(group, 13, 13)]
psa_percents_zip[ , zip3 := substr(group, 4, 6)]

dt <- rbind(
  psa_percents_zip[year == "2017" & season == "H", c("year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_zip[year == "2018" & season == "H", c("year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_zip[year == "2019" & season == "H", c("year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_zip[year == "2017" & season == "L", c("year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_zip[year == "2018" & season == "L", c("year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")],
  psa_percents_zip[year == "2019" & season == "L", c("year", "season", "a1", "a2", "a3", "a4", "zip3", "STUSPS")]
)


dt <- dcast(dt, zip3 + STUSPS + year ~ season, value.var = c("a1", "a2", "a3", "a4"))
dt[, WNVH := a3_H+a4_H]
dt[, WNVL := a3_L+a4_L]
dt[, ZIKVH := a1_H+a2_H]
dt[, ZIKVL := a1_L+a2_L]

dt <- dt_zip3[dt, on = "zip3"]




zip_to_state <- unique(dt[ , c("zip3", "STUSPS")])
# dt_states <- dt_zip3[, .(geom = st_union(geometry)),by=STUSPS]
# State outlines separated
# dt_outline_48 <- dt_states[ !(STUSPS %in% c("HI", "AK", "PR"))]
# dt_outline_HI <- dt_states[STUSPS == "HI"]
# dt_outline_AK <- dt_states[STUSPS == "AK"]
# dt_outline_PR <- dt_states[STUSPS == "PR"]
# Zip3 data tables separated
dt_48 <- dt[ !(STUSPS %in% c("HI", "AK", "PR"))]
dt_HI <- dt[STUSPS == "HI"]
dt_AK <- dt[STUSPS == "AK"] 
dt_PR <- dt[STUSPS == "PR"] 

# # dt_outline_48_mini <- dt_outline_48[ STUSPS %in% c("WA", "FL", "ME")]
# dt_48_mini <- dt_48[ STUSPS %in% c("WA", "FL", "ME")]



# 
us_outline <- st_read('data/gz_2010_us_outline_500k.shp')
us_outline <- data.table(us_outline)
# STATEFP is FIPS codes. PR is 72; AK is 02; HI is 15
# us_outline <- ms_simplify(us_outline, keep = 0.01)
us_outline_HI <- us_outline[R_STATEFP == "15" | L_STATEFP == "15" ]
us_outline_AK <- us_outline[R_STATEFP == "02" | L_STATEFP == "02" ]
us_outline_PR <- us_outline[R_STATEFP == "72" | L_STATEFP == "72" ]
us_outline_48 <- us_outline[!(R_STATEFP %in% c("15", "02", "72")) & !(L_STATEFP %in% c("15", "02", "72")) ]


#######.
# Plot the data
######

#Plotting functions ##############
makeplots <- function(dt_zip, dt_outline, fill, max = 1){
  return(
    ggplot()+
      geom_sf(data = dt_outline, aes(geometry = geometry), fill = NA)+
      geom_sf(data = dt_zip, aes(geometry = geometry, fill = fill), color = NA) + 
      facet_wrap(~year, ncol = 2) +
      theme(panel.grid=element_blank(),
            axis.ticks=element_blank(),
            axis.text=element_blank(),
            axis.title = element_blank(),
            legend.position=c(0.8, 0.23),
            panel.border=element_blank(),
            strip.background=element_blank(),
            strip.text=element_text(face="bold", hjust=0, size=10),
            panel.spacing = unit(0, "lines"))+
      scale_fill_gradientn(colours=c("#ffffff", brewer.pal(n=9, name="YlOrRd")),
                           na.value="#ffffff", name="% of iterations",
                           labels = percent, limits = c(0,max))
    
  )
  
}

annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
            inherit.aes = TRUE, params = list(grob = grob, 
                                              xmin = xmin, xmax = xmax, 
                                              ymin = ymin, ymax = ymax))
    }
  
gen_annotations <- function(row){
  return(annotation_custom2(
    grob = row$grobs[[1]],
    xmin = row$xmin[[1]],
    xmax = row$xmax[[1]],
    ymin = row$ymin[[1]],
    ymax = row$ymax[[1]],
    data = data.frame(year = row$year[[1]])))
}

##########

#
# WNVH
#


gg_HI_WNVH <- makeplots(dt_HI, us_outline_HI, fill = dt_HI$WNVH) + coord_sf(xlim = c(-161, -154), ylim = c(19, 23))
ggsave("Results/WNVH_HI_zip3.png", plot = gg_HI_WNVH, width = 6.5, height = 5)

gg_AK_WNVH <- makeplots(dt_AK, us_outline_AK, fill = dt_AK$WNVH) + coord_sf(xlim = c(-180, -129))
ggsave("Results/WNVH_AK_zip3.png", plot = gg_AK_WNVH, width = 6.5, height = 5)

gg_PR_WNVH <- makeplots(dt_PR, us_outline_PR, fill = dt_PR$WNVH)
ggsave("Results/WNVH_PR_zip3.png", plot = gg_PR_WNVH, width = 6.5, height = 5)

gg_48_WNVH <- makeplots(dt_48, us_outline_48, fill = dt_48$WNVH)
ggsave("Results/WNVH_48_zip3.png", plot = gg_48_WNVH, width = 6.5, height = 5)


HI_grobs = list(ggplotGrob(gg_HI_WNVH)[8,5],
            ggplotGrob(gg_HI_WNVH)[8,9],
            ggplotGrob(gg_HI_WNVH)[13,5])

AK_grobs = list(ggplotGrob(gg_AK_WNVH)[8,5],
             ggplotGrob(gg_AK_WNVH)[8,9],
             ggplotGrob(gg_AK_WNVH)[13,5])

PR_grobs = list(ggplotGrob(gg_PR_WNVH)[8,5],
             ggplotGrob(gg_PR_WNVH)[8,9],
             ggplotGrob(gg_PR_WNVH)[13,5])


insets_WNVH <- data.table(
  year = as.factor(rep(c(2017, 2018, 2019), 3)),
  grobs = c(HI_grobs,
            AK_grobs,
            PR_grobs),
  ymin = c(rep(c(23, 22, 23), each = 3)),
  ymax = c(rep(c(32, 34, 26), each = 3)),
  xmin = c(rep(c(-111, -125, -92), each = 3)),
  xmax = c(rep(c(-98, -108, -85), each = 3))
)


gg_WNVH <- makeplots(dt_48, us_outline_48, fill = dt_48$WNVH)+
  gen_annotations(insets_WNVH[1])+
  gen_annotations(insets_WNVH[2])+
  gen_annotations(insets_WNVH[3])+
  gen_annotations(insets_WNVH[4])+
  gen_annotations(insets_WNVH[5])+
  gen_annotations(insets_WNVH[6])+
  gen_annotations(insets_WNVH[7])+
  gen_annotations(insets_WNVH[8])+
  gen_annotations(insets_WNVH[9])

ggsave("Results/WNVH_zip3.png", plot = gg_WNVH, width = 6.5, height = 4)

#
# WNVL
#


gg_HI_WNVL <- makeplots(dt_HI, us_outline_HI, fill = dt_HI$WNVL) + coord_sf(xlim = c(-161, -154), ylim = c(19, 23))
ggsave("Results/WNVL_HI_zip3.png", plot = gg_HI_WNVL, width = 6.5, height = 5)

gg_AK_WNVL <- makeplots(dt_AK, us_outline_AK, fill = dt_AK$WNVL) + coord_sf(xlim = c(-180, -129))
ggsave("Results/WNVL_AK_zip3.png", plot = gg_AK_WNVL, width = 6.5, height = 5)

gg_PR_WNVL <- makeplots(dt_PR, us_outline_PR, fill = dt_PR$WNVL)
ggsave("Results/WNVL_PR_zip3.png", plot = gg_PR_WNVL, width = 6.5, height = 5)

gg_48_WNVL <- makeplots(dt_48, us_outline_48, fill = dt_48$WNVL)
ggsave("Results/WNVL_48_zip3.png", plot = gg_48_WNVL, width = 6.5, height = 5)


HI_grobs = list(ggplotGrob(gg_HI_WNVL)[8,5],
                ggplotGrob(gg_HI_WNVL)[8,9],
                ggplotGrob(gg_HI_WNVL)[13,5])

AK_grobs = list(ggplotGrob(gg_AK_WNVL)[8,5],
                ggplotGrob(gg_AK_WNVL)[8,9],
                ggplotGrob(gg_AK_WNVL)[13,5])

PR_grobs = list(ggplotGrob(gg_PR_WNVL)[8,5],
                ggplotGrob(gg_PR_WNVL)[8,9],
                ggplotGrob(gg_PR_WNVL)[13,5])


insets_WNVL <- data.table(
  year = as.factor(rep(c(2017, 2018, 2019), 3)),
  grobs = c(HI_grobs,
            AK_grobs,
            PR_grobs),
  ymin = c(rep(c(23, 22, 23), each = 3)),
  ymax = c(rep(c(32, 34, 26), each = 3)),
  xmin = c(rep(c(-111, -125, -92), each = 3)),
  xmax = c(rep(c(-98, -108, -85), each = 3))
)


gg_WNVL <- makeplots(dt_48, us_outline_48, fill = dt_48$WNVL)+
  gen_annotations(insets_WNVL[1])+
  gen_annotations(insets_WNVL[2])+
  gen_annotations(insets_WNVL[3])+
  gen_annotations(insets_WNVL[4])+
  gen_annotations(insets_WNVL[5])+
  gen_annotations(insets_WNVL[6])+
  gen_annotations(insets_WNVL[7])+
  gen_annotations(insets_WNVL[8])+
  gen_annotations(insets_WNVL[9])

ggsave("Results/WNVL_zip3.png", plot = gg_WNVL, width = 6.5, height = 4)

#
# ZIKVH
#


gg_HI_ZIKVH <- makeplots(dt_HI, us_outline_HI, fill = dt_HI$ZIKVH, max = 0.08) + coord_sf(xlim = c(-161, -154), ylim = c(19, 23))
ggsave("Results/ZIKVH_HI_zip3.png", plot = gg_HI_ZIKVH, width = 6.5, height = 5)

gg_AK_ZIKVH <- makeplots(dt_AK, us_outline_AK, fill = dt_AK$ZIKVH, max = 0.08) + coord_sf(xlim = c(-180, -129))
ggsave("Results/ZIKVH_AK_zip3.png", plot = gg_AK_ZIKVH, width = 6.5, height = 5)

gg_PR_ZIKVH <- makeplots(dt_PR, us_outline_PR, fill = dt_PR$ZIKVH, max = 0.08)
ggsave("Results/ZIKVH_PR_zip3.png", plot = gg_PR_ZIKVH, width = 6.5, height = 5)

gg_48_ZIKVH <- makeplots(dt_48, us_outline_48, fill = dt_48$ZIKVH, max = 0.08)
ggsave("Results/ZIKVH_48_zip3.png", plot = gg_48_ZIKVH, width = 6.5, height = 5)


HI_grobs = list(ggplotGrob(gg_HI_ZIKVH)[8,5],
                ggplotGrob(gg_HI_ZIKVH)[8,9],
                ggplotGrob(gg_HI_ZIKVH)[13,5])

AK_grobs = list(ggplotGrob(gg_AK_ZIKVH)[8,5],
                ggplotGrob(gg_AK_ZIKVH)[8,9],
                ggplotGrob(gg_AK_ZIKVH)[13,5])

PR_grobs = list(ggplotGrob(gg_PR_ZIKVH)[8,5],
                ggplotGrob(gg_PR_ZIKVH)[8,9],
                ggplotGrob(gg_PR_ZIKVH)[13,5])


insets_ZIKVH <- data.table(
  year = as.factor(rep(c(2017, 2018, 2019), 3)),
  grobs = c(HI_grobs,
            AK_grobs,
            PR_grobs),
  ymin = c(rep(c(23, 22, 23), each = 3)),
  ymax = c(rep(c(32, 34, 26), each = 3)),
  xmin = c(rep(c(-111, -125, -92), each = 3)),
  xmax = c(rep(c(-98, -108, -85), each = 3))
)


gg_ZIKVH <- makeplots(dt_48, us_outline_48, fill = dt_48$ZIKVH)+
  gen_annotations(insets_ZIKVH[1])+
  gen_annotations(insets_ZIKVH[2])+
  gen_annotations(insets_ZIKVH[3])+
  gen_annotations(insets_ZIKVH[4])+
  gen_annotations(insets_ZIKVH[5])+
  gen_annotations(insets_ZIKVH[6])+
  gen_annotations(insets_ZIKVH[7])+
  gen_annotations(insets_ZIKVH[8])+
  gen_annotations(insets_ZIKVH[9])

ggsave("Results/ZIKVH_zip3.png", plot = gg_ZIKVH, width = 6.5, height = 4)





#
# ZIKVL
#


gg_HI_ZIKVL <- makeplots(dt_HI, us_outline_HI, fill = dt_HI$ZIKVL, max = 0.25) + coord_sf(xlim = c(-161, -154), ylim = c(19, 23))
ggsave("Results/ZIKVL_HI_zip3.png", plot = gg_HI_ZIKVL, width = 6.5, height = 5)

gg_AK_ZIKVL <- makeplots(dt_AK, us_outline_AK, fill = dt_AK$ZIKVL, max = 0.25) + coord_sf(xlim = c(-180, -129))
ggsave("Results/ZIKVL_AK_zip3.png", plot = gg_AK_ZIKVL, width = 6.5, height = 5)

gg_PR_ZIKVL <- makeplots(dt_PR, us_outline_PR, fill = dt_PR$ZIKVL, max = 0.25)
ggsave("Results/ZIKVL_PR_zip3.png", plot = gg_PR_ZIKVL, width = 6.5, height = 5)

gg_48_ZIKVL <- makeplots(dt_48, us_outline_48, fill = dt_48$ZIKVL, max = 0.25)
ggsave("Results/ZIKVL_48_zip3.png", plot = gg_48_ZIKVL, width = 6.5, height = 5)


HI_grobs = list(ggplotGrob(gg_HI_ZIKVL)[8,5],
                ggplotGrob(gg_HI_ZIKVL)[8,9],
                ggplotGrob(gg_HI_ZIKVL)[13,5])

AK_grobs = list(ggplotGrob(gg_AK_ZIKVL)[8,5],
                ggplotGrob(gg_AK_ZIKVL)[8,9],
                ggplotGrob(gg_AK_ZIKVL)[13,5])

PR_grobs = list(ggplotGrob(gg_PR_ZIKVL)[8,5],
                ggplotGrob(gg_PR_ZIKVL)[8,9],
                ggplotGrob(gg_PR_ZIKVL)[13,5])


insets_ZIKVL <- data.table(
  year = as.factor(rep(c(2017, 2018, 2019), 3)),
  grobs = c(HI_grobs,
            AK_grobs,
            PR_grobs),
  ymin = c(rep(c(23, 22, 23), each = 3)),
  ymax = c(rep(c(32, 34, 26), each = 3)),
  xmin = c(rep(c(-111, -125, -92), each = 3)),
  xmax = c(rep(c(-98, -108, -85), each = 3))
)


gg_ZIKVL <- makeplots(dt_48, us_outline_48, fill = dt_48$ZIKVL)+
  gen_annotations(insets_ZIKVL[1])+
  gen_annotations(insets_ZIKVL[2])+
  gen_annotations(insets_ZIKVL[3])+
  gen_annotations(insets_ZIKVL[4])+
  gen_annotations(insets_ZIKVL[5])+
  gen_annotations(insets_ZIKVL[6])+
  gen_annotations(insets_ZIKVL[7])+
  gen_annotations(insets_ZIKVL[8])+
  gen_annotations(insets_ZIKVL[9])

ggsave("Results/ZIKVL_zip3.png", plot = gg_ZIKVL, width = 6.5, height = 4)






# 
# 
# 
# usa_zip5 <- usa_zip5[ , c("ZCTA5CE10", "geometry")]
# usa_zip5[ , ZTCA3 := substr(ZCTA5CE10, 1, 3)]
# 
# ggplot()+
#   geom_sf(data = usa_zip3[1:20, ], aes(geometry = geometry))
# 


# psa_percents <- fread("Results/PSA_opt_by_state_combined.csv")
# basecase_results[ , year := substr(group, start=1, stop = 4)]
# basecase_results[ , STUSPS := substr(group, start=6, stop = 7)]
# basecase_results$group
#
#
#
# usa_48 <- usa[!(usa$STUSPS %in% c("PR", "AK", "HI")) , ]$STUSPS
#
# usa_48_plt <- ggplot(data=usa[usa$STUSPS %in% usa_48, ])+
#   geom_sf()+
#   theme_void()+
#   ylim(20, 50)
#
# #usa_48_plt
#
# usa_AK <- ggplot(data=usa[usa$STUSPS == "AK",])+
#   geom_sf()+
#   theme_void()+
#   coord_sf(xlim = c(-180, -129))
#
# #usa_AK
#
# usa_HI <- ggplot(data=usa[usa$STUSPS == "HI",])+
#   geom_sf()+
#   coord_sf(xlim = c(-161, -154))+
#   theme_void()
#
# #usa_HI
#
# usa_PR <- ggplot(data=usa[usa$STUSPS == "PR",])+
#   geom_sf()+
#   theme_void()
#
# grob_2017 <- ggplotGrob(ggdraw(usa_48_plt)+
#                           draw_plot(usa_AK, x=0.03, y=0.16, width = 0.2, height = 0.2)+
#                           draw_plot(usa_HI, x=0.22, y=0.16, width = 0.2, height = 0.2)+
#                           draw_plot(usa_PR, x=0.55, y=0.18, width = 0.1, height = 0.1))
#
#
# # Read in results
# basecase_results <- fread("Results/opt_by_state_basecase.csv")
# psa_metrics <- fread("Results/opt_by_iter_metric_allyears_combined.csv")
# psa_percents <- fread("Results/opt_by_group_combined.csv")
#
#
#
# # 2019
#
# labs2019_h<- data.table(statepop)
# labs2019_h[ , pop_2015 := NULL]
# labs2019_h[ , Policy := as.factor(c("No intervention",
#                               "No intervention",
#                               "WNV MP-NAT",
#                               rep("No intervention", 51-3)))]
#
# p_2019_h <- plot_usmap(regions = "states",
#            data = labs2019_h,
#            values = "Policy",
#            color = "white")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c("#999999", "#8D1516", "#56B4E9"),
#                     name = "")
#
#
# labs2019_l<- data.table(statepop)
# labs2019_l[ , pop_2015 := NULL]
# labs2019_l[ , Policy := as.factor(c(rep("No intervention", 51)))]
#
# p_2019_l <- plot_usmap(regions = "states",
#                        data = labs2019_l,
#                        values = "Policy",
#                        color = "white")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c("#999999", "#8D1516", "#56B4E9"),
#                     name = "")
#
#
# # ggarrange(p_2019_h, p_2019_l,
# #           ncol = 2, common.legend = TRUE, legend="bottom",
# #           labels = c("High misquito season", "Low mosquito season"))
#
# #  2018
#
# labs2018_h<- data.table(statepop)
# labs2018_h[ , pop_2015 := NULL]
# labs2018_h[ , Policy := as.factor(c(rep("No intervention", 8),
#                                     "WNV MP-NAT", #DC
#                                     rep("No intervention", 6),
#                                     "WNV MP-NAT", #IA
#                                     rep("No intervention", 2),
#                                     "WNV MP-NAT", #LA
#                                     rep("No intervention", 7),
#                                     "WNV MP-NAT", "WNV MP-NAT", #MT NE
#                                     rep("No intervention", 6),
#                                     "WNV MP-NAT", #ND
#                                     rep("No intervention", 6),
#                                     "WNV MP-NAT", #SD
#                                     rep("No intervention", 9)))]
#
#
# dc_data <- data.frame(lon = -77.03637, lat = 38.89511)
# transformed_data <- usmap_transform(dc_data)
#
#
# p_2018_h <- plot_usmap(regions = "states",
#                        data = labs2018_h,
#                        values = "Policy",
#                        color = "white")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c("#999999", "#8D1516", "#56B4E9"),
#                     name = "") +
#   geom_point(data = transformed_data,
#              aes(x = lon.1, y = lat.1),
#              color = "#8D1516",
#              size = 3)
#
#
# labs2018_l<- data.table(statepop)
# labs2018_l[ , pop_2015 := NULL]
# labs2018_l[ , Policy := as.factor(c(rep("No intervention", 51)))]
#
# p_2018_l <- plot_usmap(regions = "states",
#                        data = labs2018_l,
#                        values = "Policy",
#                        color = "white")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c("#999999", "#8D1516", "#56B4E9"),
#                     name = "")
#
#
#
#
# # 2017
#
#
# labs2017_h<- data.table(statepop)
# labs2017_h[ , pop_2015 := NULL]
# labs2017_h[ , Policy := as.factor(c(rep("No intervention", 24),
#                                     "WNV MP-NAT", #MS
#                                     rep("No intervention", 2),
#                                     "WNV MP-NAT", "WNV MP-NAT", #NE NV
#                                     rep("No intervention", 5),
#                                     "WNV MP-NAT", #ND
#                                     rep("No intervention", 6),
#                                     "WNV MP-NAT", #SD
#                                     rep("No intervention", 2),
#                                     "WNV MP-NAT", #UT
#                                     rep("No intervention", 6)))]
#
# p_2017_h <- plot_usmap(regions = "states",
#                        data = labs2017_h,
#                        values = "Policy",
#                        color = "white")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c("#999999", "#8D1516", "#56B4E9"),
#                     name = "")
#
#
# labs2017_l<- data.table(statepop)
# labs2017_l[ , pop_2015 := NULL]
# labs2017_l[ , Policy := as.factor(c(rep("No intervention", 51)))]
#
# p_2017_l <- plot_usmap(regions = "states",
#                        data = labs2017_l,
#                        values = "Policy",
#                        color = "white")+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=c("#999999", "#8D1516", "#56B4E9"),
#                     name = "")
#
#
# fig <- ggarrange(p_2017_h, p_2018_h, p_2019_h, p_2017_l, p_2018_l, p_2019_l,
#           ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom",
#           labels = c("2017", "2018", "2019"))
#
# fig_annotated <- annotate_figure(fig, left = text_grob("          Low mosquito                          High mosquito", face = "bold", rot = 90, size = 14))
#
# ggsave("policy map.jpg", width = 10, height = 6, units = "in")
















# 
# 
# 
# 
# # Format map
# us_map <- geojson_read("us_states_hexgrid.geojson",  what = "sp")
# us_map <- data.table(fortify(us_map, region="iso3166_2"))
# 
# # Add in PR
# PR <- data.table(us_map[us_map$id=="FL", ])
# PR <- PR[, id := "PR"]
# PR <- PR[, group := "PR.1"]
# PR <- PR[, order := (max(us_map$order)+1):(max(us_map$order) + nrow(PR))]
# PR <- PR[, long := long + 12]
# 
# us_map <- rbind(us_map, PR)
# 
# # Get centers
# centers <- us_map[ , list(x = mean(long), y = mean(lat)), by = id]
# 
# # Data for map
# psa_percents <- fread("Results/PSA_opt_by_state_combined.csv")
# basecase_results[ , year := substr(group, start=1, stop = 4)]
# basecase_results[ , STUSPS := substr(group, start=6, stop = 7)]
# 
# psa_percents <- fread("Results/PSA_opt_by_group_combined.csv")[ , lapply(.SD, sum), .SD = 3:8, by = c("group", "year")]
# psa_percents[ , z := (10000 - z)/10000]
# psa_percents[, c("a1", "a2", "a3", "a4", "m1") :=lapply(.SD, function(x) x/10000), .SDcols=4:8]
# psa_percents[ , STUSPS := substr(group, 1, 2)]
# psa_percents[ , season := substr(group, 4, 4)]
# WNV_H <- psa_percents[season == "H" , value := (a3+a4)]
# 
