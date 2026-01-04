# SVG files are from:
# https://www.phylopic.org/images/938008ab-0501-4c4d-88c5-fee17189ead3/bos-primigenius-taurus
# https://www.phylopic.org/images/008d6d88-d1be-470a-8c70-73625c3fb4fb/sus-scrofa-domestica
# https://www.phylopic.org/images/3505388e-32ed-464d-816a-0bfca1b68934/meleagris-gallopavo
# https://www.phylopic.org/images/f16a316c-6580-4cc0-b6b5-fd2e895ba225/gallus-gallus-domesticus

require(sf)
require(ggplot2)
require(tigris)
require(fields)
require(dplyr)
require(gridExtra)
library(grid)
library(cowplot)
library(magick)

nvariables <- 12
input_data <- read.csv('../Datasets2024/April11_2025/seasonal_amr_noGuam.csv')

# Create a subset removing other variables not listed in table 1 of the manuscript:
table_1_Ncounts <- aggregate(.~establishment_number,data=input_data %>% select(-c(id, season, state, Sum_resistants, establishment_name, grant_date, activities, dbas, latitude, longitude, county)),FUN=sum)
for (i in 1:ncol(table_1_Ncounts)) 
    {
    table_1_Ncounts[,i] <- ifelse(table_1_Ncounts[,i]>0,1,0)
    }

ncounts <- 1:ncol(table_1_Ncounts)
for (i in 1:ncol(table_1_Ncounts))
    {
    ncounts[i] <- sum(table_1_Ncounts[,i])
    }
names(ncounts) <- colnames(table_1_Ncounts)

gnnwr_results <- read.csv(file='../gnnwr/demo/seasonal_amr_noGuam_result.csv', header=TRUE)
locations <- unique(input_data[,c('latitude','longitude')])

# Determine the weighted contribution of a particular variable for that record id (i.e., beta_j * xij, where xij is an indicator variable):

coefficients_matrix <- data.frame(matrix(nrow=nrow(gnnwr_results), ncol=nvariables))
colnames(coefficients_matrix) <- colnames(gnnwr_results)[1:nvariables]
reordered_results_matrix <- data.frame(matrix(nrow=nrow(gnnwr_results), ncol=2)) # predicted and residual
indicator_matrix <- coefficients_matrix

for (i in 1:nvariables)
    {
    indicator_matrix[,i] <- ifelse(input_data[,12+i]==0,0,1)
    }

for (i in 1:nrow(input_data))
    {
    results_row <- which(gnnwr_results[,'id']==input_data[i,'id'])
    coefficients_matrix[i,1:nvariables] <- c(indicator_matrix[i,] * gnnwr_results[results_row,1:nvariables])
    reordered_results_matrix[i,1] <- log(gnnwr_results[results_row,'Pred_Sum_resistants'])
    reordered_results_matrix[i,2] <- log(abs(gnnwr_results[results_row,'Pred_Sum_resistants'] - input_data[i,'Sum_resistants']))
    if (i %% 2500==0)
        {
        print(paste('rearranging datapoint',i, 'with id', input_data[i,'id'],'completed'))
        }
    }

coefs_by_id <- cbind(reordered_results_matrix,coefficients_matrix, input_data[,c('season','longitude','latitude')])
plot_titles <- c('Log Predicted AbR', 'Log Abs Residual AbR','Plant Size: Large', 'Plant Size: Small', 'Plant Size: Very Small', 'Chicken Broiler / Young Carcass', gsub('  ',' ', gsub('Comminuted or Otherwise Nonintact','Com./O.W.N.', gsub('_',' ',gsub('coef_Product_','',colnames(gnnwr_results)[5:(nvariables-1)])))), 'Turkey Carcass Sponge')

# Stolen from https://stackoverflow.com/a/65233844/1790399
my_sf <- st_as_sf(coefs_by_id, coords = c('longitude', 'latitude'))
st_crs(my_sf) <- 4326 # use ESPG 4326, i.e., WGS84 coordinate system

options(tigris_use_cache = TRUE, tigris_class = "sf")

us <- states(cb = TRUE) %>% filter(STUSPS != "GU" & STUSPS != "AS" & STUSPS != 'VI' & STUSPS != 'MP') %>%
    shift_geometry() %>%
    st_transform(4326) # Transform to same CRS as my_sf

draw_map <- function(coefficient, season, panel_label)
    {
    sez <- which(my_sf$season==season)
    column_in_input <- 12+coefficient-2 # the -2 is needed to account for the predicted values and the residuals in my_sf
    if(coefficient > 2)
        {
        ignoreds <- which(input_data[sez,column_in_input]==0) 
        filtered_coefficients <- my_sf[sez[-ignoreds],]
        }
    else
        {
        filtered_coefficients <- my_sf[sez,]
        }

    filtered_coefficients <- shift_geometry(filtered_coefficients)

    title_text <- paste(LETTERS[panel_label], plot_titles[coefficient], 'in', season)
    # Break after 40 characters or at a specific point
    if (nchar(title_text) > 40) {
    # Find a space near the midpoint to break at
    midpoint <- nchar(title_text) / 2
    space_positions <- gregexpr(" ", title_text)[[1]]
    best_break <- space_positions[which.min(abs(space_positions - midpoint))]
    title_text <- paste0(
        substr(title_text, 1, best_break-1), 
        "\n",
        substr(title_text, best_break+1, nchar(title_text))
    )
    }
    labs(title = title_text)

    ggplot(filtered_coefficients) + 
        geom_sf(data = us, fill = "white", color = "gray", linewidth = 0.2) +
        # overlay points on states
        geom_sf(data=filtered_coefficients[,coefficient]) + 
        # overlay color scheme to be the range of 0the coefficient effects; the 2+ is needed to account for the predictions and residuals
        geom_sf(aes(color = .data[[colnames(my_sf)[coefficient]]]),inherit.aes=FALSE) +
        scale_color_gradientn(colors = tim.colors(100), limits = range(my_sf[[colnames(my_sf)[coefficient]]][-which(input_data[,]==0)], na.rm = TRUE) ) +
        labs(title = title_text, color = ifelse(coefficient>2,'Coefficient','Value')) +
        theme(
            plot.title = element_text(
                size = 11,               # Adjust title size (default is usually 14-16)
                hjust = 0.5,            # Center title (0 = left, 0.5 = center, 1 = right)
                margin = margin(b = 10) # Add bottom margin below title
            ),
            plot.margin = margin(t = 5, r = 10, b = 5, l = 10) # Adjust overall plot margins
        )
    }

seasons_to_draw <- unique(coefs_by_id[,'season'])[seq(1,13,by=3)] # use five interspersed records

# model performance : predicted and residual
plots <- list()

k <- 1
for (i in 1:2)
    {
    for (j in 1:length(seasons_to_draw))
        {
        plots[[k]] <-draw_map(i, seasons_to_draw[j], k)
        k <- k + 1
        }
    }

pdf(file='Fig2.pdf', width = 18.958333, height = 5.277)
grid_plot <- grid.arrange(
    grobs = plots,
    ncol = 5,
    nrow = 2,
    top = textGrob("Figure 2", 
                 gp = gpar(fontsize = 16, fontface = "bold"),
                 just = "center",
                 vjust = 1)
)
dev.off()

# Impact of size
plots <- list()
k <- 1
for (i in 3:5)
    {
    for (j in 1:length(seasons_to_draw))
        {
        plots[[k]] <-draw_map(i, seasons_to_draw[j], k)
        k <- k + 1
        }
    }
pdf(file='Fig3.pdf', width = 19.229167, height = 6.555665)
grid_plot <- grid.arrange(
    grobs = plots,
    ncol = 5,
    nrow = 3,
    top = textGrob("Figure 3", 
                gp = gpar(fontsize = 16, fontface = "bold"),
                just = "center",
                vjust = 1)
)
dev.off()

# Impact of chickens
plots <- list()
k <- 1
for (i in 6:8)
    {
    for (j in 1:length(seasons_to_draw))
        {
        plots[[k]] <-draw_map(i, seasons_to_draw[j], k)
        k <- k + 1
        }
    }

grid_plot <- arrangeGrob(
  grobs = plots,
  ncol = 5,
  nrow = 3,
    top = textGrob("Figure 4", 
                gp = gpar(fontsize = 16, fontface = "bold"),
                just = "center",
                vjust = 1)
)

# Get the taxon silhoutte. Note to use the silhoutte with svg
# you need to install imagemagick with rsvg or else it looks bad.
# via sudo apt-get install librsvg2-dev
# Then compile imagemagick from source with rsvg support
# https://imagemagick.org/script/install-source.php#gsc.tab=0 via
# git clone --depth 1 https://github.com/ImageMagick/ImageMagick.git ImageMagick-7.1.2
# ./configure --with-rsvg=yes --enable-shared
# make
# sudo make install

draw_with_svg <- ggdraw() +
    draw_plot(grid_plot) +
    draw_image('chicken.svg', x=0.9, y=0.89, width=0.095, height=0.095)

pdf(file='Fig4_with_svg.pdf', width = 19.229167, height = 6.555665)
draw_with_svg
dev.off()


# Impact of beef and pork
plots <- list()
k <- 1
for (i in 9:12)
    {
    for (j in 1:length(seasons_to_draw))
        {
        plots[[k]] <-draw_map(i, seasons_to_draw[j], k)
        k <- k + 1
        }
    }

grid_plot <- arrangeGrob(
  grobs = plots,
  ncol = 5,
  nrow = 4,
  padding = unit(0.1, "line")
)

# Wrap in cowplot to add top margin
grid_plot <- ggdraw() +
  draw_plot(grid_plot, y = 0, height = 0.90)  # Leave 10% space at top

draw_with_svg <- grid_plot +
    draw_label("Figure 5", x = 0.5, y = 0.95, size = 18, fontface = "bold") +
    draw_image('cow.svg', x=0.8, y=0.89, width=0.095, height=0.095) + 
    draw_label("/", x=0.895, y=0.93, size=35, fontface="bold") +
    draw_image('swine.svg', x=0.9, y=0.89, width=0.08, height=0.08)

pdf(file='Fig5_with_svg.pdf', width = 19.229167, height = 9.215332)
draw_with_svg
dev.off()

# Turkey
plots <- list()
k <- 1
for (i in 13:14)
    {
    for (j in 1:length(seasons_to_draw))
        {
        plots[[k]] <-draw_map(i, seasons_to_draw[j], k)
        k <- k + 1
        }
    }
grid_plot <- arrangeGrob(
  grobs = plots,
  ncol = 5,
  nrow = 2,
    top = textGrob("Figure 6", 
                gp = gpar(fontsize = 16, fontface = "bold"),
                just = "center",
                vjust = 1)
)

draw_with_svg <- ggdraw() +
    draw_plot(grid_plot) +
    draw_image('turkey.svg', x=0.9, y=0.89, width=0.095, height=0.095)

pdf(file='Fig6_with_svg.pdf', width = 19.229167, height = 6.555665)
draw_with_svg
dev.off()

# Supplementary figures for all seasons
for (coefficient in 1:(nvariables+2)) # +2 to account for the predicted values and the residuals
    {
    plots <- list()
    k <- 1
    for (season in unique(coefs_by_id[,'season']))
        {
        plots[[k]] <- draw_map(coefficient, season, k)
        k <- k + 1
        }

    pdf(file=paste('S1',letters[coefficient],'.pdf',sep=''), width = 19.229167, height = 6.555665)
    grid_plot <- grid.arrange(
    grobs = plots,
    layout_matrix = rbind(
        matrix(1:5, ncol=5),
        matrix(6:10, ncol=5),
        cbind(matrix(NA, ncol=1), matrix(11:13, ncol=3), matrix(NA, ncol=1))
    ),
        top = textGrob(plot_titles[coefficient], 
                    gp = gpar(fontsize = 16, fontface = "bold"),
                    just = "center",
                    vjust = 1)
    )
    dev.off()
    }

pdf(file='S2.pdf')
run_stats <- read.delim(file='../gnnwr/gtnnwr_logs/gtnnwr20250731-232913.log',sep=' ',header=F)
validLoss <- as.numeric( gsub(';','',run_stats[,20]))
plot(validLoss,t='l',log='x',xlab='Epoch',ylab='Validation Loss')
dev.off()
