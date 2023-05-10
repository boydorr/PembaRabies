shp_to_df <- function(shapefile){

  # Transform shapefile for plotting
  shapefile@data$id <- rownames(shapefile@data)
  shp_df <- fortify(shapefile, region = "id")
  shp_df <- merge(shp_df, shapefile@data, by = "id")

  return(shp_df)
}
