get_FCS_metadata <- function(paths){
  # Function to retrieve metadata parameters from a list of FCS paths.
  # Args:
  #   paths: list of file locations.
  # Returns:
  #   dataframe containing metadata information on every file used for downstream
  #   QC and plotting.
  headers <- read.FCSheader(paths)
  # Iterate over the headers and convert them to dataframes
  parameter_list <- list()
  for (i in 1:length(headers)) {
    header <- headers[[i]]
    # First, identify all the markers that were used in this FCS
    # If markers were not measured, they are often excluded as $P_S parameter
    # but still have a $P_N and $P_V parameter
    p_markers <- header[grepl("^\\$P[1-9][0-9]?[S]$", names(header))]
    p_names <- header[grepl("^\\$P[1-9][0-9]?[N]$", names(header))]
    p_voltages <- header[grepl("^\\$P[1-9][0-9]?[V]$", names(header))]
    # Keep only the values for markers
    p_names <- p_names[as.character(lapply(names(p_markers), function(x) gsub("S", "N", x)))]
    p_voltages <- p_voltages[as.character(lapply(names(p_markers), function(x) gsub("S", "V", x)))]
    
    # Add all the variables to a dataframe
    values <- list()
    values['file'] <- header[['$FIL']]
    values['GUID'] <- header[['GUID']]
    values['time'] <- header[['EXPORT TIME']]
    values['cytometer_id'] <- header[['CYTNUM']]
    values['channel_str'] <- paste(p_names, collapse = "_")
    values['marker_str'] <- paste(p_markers, collapse = "_")
    values['voltage_str'] <- paste(p_voltages, collapse = "_")
    parameter_list[[i]] <- values
  }
  parameters <- dplyr::bind_rows(parameter_list)
  # Create a new shorter categorical representation of the voltage string
  parameters <- parameters[order(parameters$voltage_str),]
  parameters$voltage_group <- as.integer(factor(parameters$voltage_str))
  parameters$voltage_group <- paste("VOLTAGEGROUP", parameters$voltage_group, sep = "")
  # Create a new shorter categorical representation of the cytometers
  parameters <- parameters[order(parameters$cytometer_id),]
  parameters$cytometer_group <- as.integer(factor(parameters$cytometer_id))
  parameters$cytometer_group <- paste("CYTOMETER", parameters$cytometer_group, sep = "")
  # Create a combined representation
  parameters$cytometer_voltage <- apply(parameters[,c("cytometer_group", "voltage_group")], 1, function(x) paste(x, collapse = "_"))
  parameters <- data.frame(parameters)
  return(parameters)
}