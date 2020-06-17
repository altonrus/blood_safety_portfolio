library(ggplot2)
library(data.table)
library(mapproj)
library(scales)
#library(sf)
#library(cowplot)
#library(gridExtra)
#library(rgeos)
theme_set(theme_minimal())

library(geojsonio)


# Format map
us_map <- geojson_read("us_states_hexgrid.geojson",  what = "sp")
us_map <- data.table(fortify(us_map, region="iso3166_2"))

# Add in PR
PR <- data.table(us_map[us_map$id=="FL", ])
PR <- PR[, id := "PR"]
PR <- PR[, group := "PR.1"]
PR <- PR[, order := (max(us_map$order)+1):(max(us_map$order) + nrow(PR))]
PR <- PR[, long := long + 12]

us_map <- rbind(us_map, PR)

# Get centers
centers <- us_map[ , list(x = mean(long), y = mean(lat)), by = id]

# Data for map
psa_percents <- fread("Results/PSA_opt_by_state_combined.csv")
basecase_results[ , year := substr(group, start=1, stop = 4)]
basecase_results[ , STUSPS := substr(group, start=6, stop = 7)]

psa_percents <- fread("Results/PSA_opt_by_group_combined.csv")[ , lapply(.SD, sum), .SD = 3:8, by = c("group", "year")]
psa_percents[ , z := (10000 - z)/10000]
psa_percents[, c("a1", "a2", "a3", "a4", "m1") :=lapply(.SD, function(x) x/10000), .SDcols=4:8]
psa_percents[ , STUSPS := substr(group, 1, 2)]
psa_percents[ , season := substr(group, 4, 4)]
WNV_H <- psa_percents[season == "H" , value := (a3+a4)]



# Gen fig

gg <- ggplot()
gg <- gg + geom_map(data=us_map, map=us_map,
                    aes(x=long, y=lat, map_id=id),
                    color="white", size=0.5)
gg <- gg + geom_map(data=WNV_H, map=us_map,
                    aes(fill=value, map_id=STUSPS))
gg <- gg + geom_map(data=WNV_H, map=us_map,
                    aes(fill=value, map_id=STUSPS), alpha = 0, color = "white")
gg <- gg + geom_text(data=centers, aes(label=id, x=x, y=y), color="white", size=3)
gg <- gg + coord_map()
gg <- gg + facet_wrap(~year, ncol = 2)
gg <- gg + theme(panel.grid=element_blank())
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text=element_blank())
#gg <- gg + theme(legend.direction="horizontal")
gg <- gg + theme(axis.title = element_blank())
#gg <- gg + labs(fill = "% iterations WNV testing optimal")
gg <- gg + theme(legend.position=c(0.8, 0.23))
gg <- gg + theme(panel.border=element_blank())
gg <- gg + guides(fill = guide_legend())
gg <- gg + theme(strip.background=element_blank())
gg <- gg + theme(strip.text=element_text(face="bold", hjust=0, size=10))
gg <- gg + theme(panel.spacing = unit(0, "lines"))
gg <- gg + scale_fill_gradient(name = "", #name="% iterations WNV\ntesting optimal",
                               high = "#132B43",
                               low = "#56B1F7",
                               labels = percent
                               )

ggsave("Results/WNV_policy_hexbin.png", plot = gg, width = 6.5, height = 5)
gg <- gg + ggtitle("Percent of probabilistic sensitivity analysis iterations for which West Nile Virus
                   testing during high mosquito season was optimal by area")
ggsave("Results/WNV_policy_hexbin_wtitle.png", plot = gg, width = 6.5, height = 5)






# ggsave("SMDM_fig.png",
#        plot = grid.arrange(
#              p_policy + ggtitle("A. Optimal blood safety portfolio as a function of Zika and West Nile virus prevalence"),
#              gg + ggtitle("B. Percent of probabilistic sensitivity analysis iterations for which West Nile Virus testing during high mosquito season was optimal by state"),#+
#                #theme(plot.title = element_text(hjust = 1)), 
#              ncol = 1,
#              widths = unit(c(6), "in"),
#              heights = unit(c(4.5, 4.5), "in")),
#        width = 6,
#        height = 9,
#        units = "in"
# )





# 
# 
# 
# 
# 
# 
# 
# usa <- st_read('cb_2018_us_state_20m/cb_2018_us_state_20m.shp')
# usa <- usa[!(usa$STUSPS %in% c("VI", "MP", "GU", "AS")) , ]
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
