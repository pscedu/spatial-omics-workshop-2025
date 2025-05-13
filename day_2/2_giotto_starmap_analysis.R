library(Giotto)

##%######################################################%##
#                                                          #
####                  Download dataset                  ####
#                                                          #
##%######################################################%##

GiottoData::getSpatialDataset(
    dataset = 'starmap_3D_cortex',
    directory = "./",
    method = 'wget'
)

##%######################################################%##
#                                                          #
####                Create Giotto object                ####
#                                                          #
##%######################################################%##

instructions <- createGiottoInstructions(save_plot = FALSE,
                                         show_plot = TRUE,
                                         return_plot = FALSE,
                                         python_path = NULL)

g <- createGiottoObject(expression = "STARmap_3D_data_expression.txt",
                        spatial_locs = "STARmap_3D_data_cell_locations.txt",
                        instructions = instructions)

##%######################################################%##
#                                                          #
####                Visualize 3D object                 ####
#                                                          #
##%######################################################%##

spatPlot3D(g,
           point_size = 1)

##%######################################################%##
#                                                          #
####                     Filtering                      ####
#                                                          #
##%######################################################%##

g <- filterGiotto(g,
                  expression_threshold = 1,
                  feat_det_in_min_cells = 1,
                  min_det_feats_per_cell = 10)

##%######################################################%##
#                                                          #
####                   Normalization                    ####
#                                                          #
##%######################################################%##

g <- normalizeGiotto(g)

##%######################################################%##
#                                                          #
####                     Statistics                     ####
#                                                          #
##%######################################################%##

g <- addStatistics(g)

spatPlot3D(g,
           point_size = 1,
           cell_color = "total_expr",
           color_as_factor = FALSE)

# By default, this is displayed in 3D "cube" scaling. 
# We can change this to "real" to get an undistorted view.
spatPlot3D(g,
           point_size = 1,
           cell_color = "total_expr",
           color_as_factor = FALSE,
           axis_scale = "real")

##%######################################################%##
#                                                          #
####                Dimension reduction                 ####
#                                                          #
##%######################################################%##

g <- runPCA(g)

##%######################################################%##
#                                                          #
####                        UMAP                        ####
#                                                          #
##%######################################################%##

g <- runUMAP(g,
             dimensions_to_use = 1:15,
             n_components = 3 # make 3D UMAP
)

##%######################################################%##
#                                                          #
####                         NN                         ####
#                                                          #
##%######################################################%##

g <- createNearestNetwork(g, 
                          dimensions = 1:15)

##%######################################################%##
#                                                          #
####                 Leiden Clustering                  ####
#                                                          #
##%######################################################%##

g <- doLeidenCluster(g, 
                     resolution = 0.2)

##%######################################################%##
#                                                          #
####                      3D plots                      ####
#                                                          #
##%######################################################%##

plotUMAP_3D(g,
            point_size = 1,
            cell_color = "leiden_clus", 
            show_center_label = FALSE)

spatPlot3D(g,
           point_size = 1,
           cell_color = "leiden_clus")

spatPlot3D(g,
           point_size = 1,
           cell_color = "leiden_clus",
           axis_scale = "real")




