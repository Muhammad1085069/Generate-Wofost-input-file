
generate_wofost_rotations_multi <- function(rotations,lon_lat_pairs) {
# Function to adjust years for each rotation
adjust_years <- function(rotation) {
  years <- rep(0, length(rotation))
  for (i in 1:length(rotation)) {
    if (i == 1) {
      years[i] <- 2020
    } else {
      if (rotation[i] == "fallow" && rotation[i - 1] %in% c("wheat", "rapeseed")) {
        years[i] <- years[i - 1] + 1  # Fallow in the next year after wheat or rapeseed
      } else if (rotation[i] == "fallow" && rotation[i - 1] %in% c("maize", "barley", "potato", "rye")) {
        years[i] <- years[i - 1]  # Fallow in the same year as starting year for maize, barley, potato, or rye
      } else if (rotation[i] == "fallow") {
        years[i] <- years[i - 1] + 1  # Fallow in the next year
      } else if (rotation[i] %in% c("barley", "maize", "potato", "rye")) {
        years[i] <- years[i - 1] + 1  # Crops after fallow in the same year
      } else {
        years[i] <- years[i - 1]  # Default: Same year
      }
    }
  }
  return(years)
}

# Apply the adjust_years function to each rotation
rotation_years <- lapply(rotations, adjust_years)
# Initialize a list to store rotation data frames
rotation_dfs <- list()

# Iterate over rotations and lon-lat pairs
for (rotation_name in names(rotations)) {
  for (i in 1:nrow(lon_lat_pairs)) {
    lon <- lon_lat_pairs$lon[i]
    lat <- lon_lat_pairs$lat[i]
    iteration <- i
    rotation <- rotations[[rotation_name]]
    years <- rotation_years[[rotation_name]]
    
    # Create a rotation dataframe for each lon-lat pair
    rotation_df <- data.frame(
      #rotation_id = paste0(rotation_name, "_", i, ".", seq_along(rotation)),
      rotation = rep(paste0(rotation_name, ".", i), length(rotation)),
      crop = rotation,
      year = years,
      lon = rep(lon, length(rotation)),
      lat = rep(lat, length(rotation)),
      iteration = rep(iteration, length(rotations[[rotation_name]]))
    )
    
    rotation_dfs[[paste0(rotation_name, "_", i)]] <- rotation_df
  }
}

combined_df <- do.call(rbind, rotation_dfs)
combined_df <- combined_df[order(combined_df$iteration), ]

# Define a function to add suffix to crop names
add_variety <- function(crop) {
  if (crop == "wheat") {
    return("Winter_wheat_101")
  } else if (crop == "maize") {
    return("Grain_maize_201")
  } else if (crop == "barley") {
    return("Spring_barley_301")
  } else if (crop == "rye_grass") {
    return("Northern_RyeGrass")
  } else if (crop == "rapeseed") {
    return("Oilseed_rape_1001")
  } else if (crop == "potato") {
    return("Potato_701")
  } else if (crop == "fallow") {
    return("fallow")
  } else {
    return(NA)
  }
}

# Apply the function to create a new column with variety names
combined_df$variety <- sapply(combined_df$crop, add_variety)
# Define column names for the LHS inputs
cols <- c("WAV", "SMLIM", "NAVAILI", "PAVAILI", "KAVAILI",
          "N_1", "N_2", "N_3", "N_4", "P_1", "P_2", "P_3","P_4", "K_1", "K_2", "K_3","K_4")

# Define minimum and maximum values for each input variable
min_values <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0)
max_values <- c(60, 0.5, 100, 50, 100, 250, 250, 250, 250, 250, 250, 250, 250, 250,250,250,250)

# Create a function to descale the inputs
descale_inputs <- function(input_n, min_values, max_values) {
  if (!is.matrix(input_n)){
    input_n<- t(input_n)
  }
  num_rows <- ifelse(is.null(nrow(input_n)), 1, nrow(input_n)) 
  num_cols <- ifelse(is.null(ncol(input_n)), length(input_n), ncol(input_n)) 
  descaled_inputs <- matrix(0, nrow = num_rows, ncol = num_cols)
  
  for (i in 1:num_rows) {
    for (j in 1:num_cols) {
      descaled_inputs[i, j] <- (input_n[i, j] * (max_values[j] - min_values[j])) + min_values[j]
    }
  }
  
  return(descaled_inputs)
}
library(lhs)
# Define the number of sets to generate
num_sets <- 1000

# Create an empty list to store the LHS samples
lhs_samplesa <- list()
descaled_samples <- list()
# Loop to generate LHS samples
for (i in 1:num_sets) {
  #df3 <- maximinLHS(nrow(data_combined), length(cols))
  df3=maximinLHS(nrow(combined_df), length(cols))
  lhs_samplesa[[i]] <- df3
  descaled_samples[[i]] <- descale_inputs(lhs_samplesa[[i]], min_values, max_values)
}

#df3 <- maximinLHS(nrow(data_combined), 17) 
input_data <- do.call(rbind, descaled_samples)
colnames(input_data) <- cols

# Function to replace empty cells for N_1 to K_4 for the crop "fallow"
replace_empty_cells <- function(df) {
  # Iterate through each row of the dataframe
  for (i in 1:nrow(df)) {
    # Check if the crop is "fallow"
    if (df[i, "crop"] == "fallow") {
      # For "fallow", clear all N_1 to K_4
      for (j in 1:4) {
        df[i, paste0("N_", j)] <- ""
        df[i, paste0("P_", j)] <- ""
        df[i, paste0("K_", j)] <- ""
      }
    } else if (df[i, "crop"] %in% c("wheat", "barley", "maize")) {
      # For "wheat", "barley", and "maize", only clear N_4, P_4, and K_4
      df[i, "N_4"] <- ""
      df[i, "P_4"] <- ""
      df[i, "K_4"] <- ""
    }
  }
  return(df)
}

# Combine LHS inputs with the combined dataframe
replicated_matrices=lapply(1:num_sets, function(x) combined_df)
data_combined <- do.call(rbind, replicated_matrices)
data_combined$iteration<- rep(1:num_sets, each =  length(combined_df$iteration ))
data_f <- cbind(input_data, data_combined)



# Generate LHS samples
set.seed(123456)
lhs_samples <- t(maximinLHS(1, nrow(data_f)))
lhs_samples=cbind(lhs_samples, data_f$year, data_f$crop)


calculate_crop_start_date <- function(year, crop) {
  if (crop == "wheat") {
    return(as.numeric(as.POSIXct(paste(year, "-11-15", sep = ""), format = "%Y-%m-%d")))
  } else if (crop == "fallow") {
    return(as.numeric(as.POSIXct(paste(year, "-09-01", sep = ""), format = "%Y-%m-%d")))
  } else if (crop == "maize") {
    return(as.numeric(as.POSIXct(paste(year, "-03-01", sep = ""), format = "%Y-%m-%d")))
  } else if (crop == "barley") {
    return(as.numeric(as.POSIXct(paste(year, "-02-15", sep = ""), format = "%Y-%m-%d")))
  } else if (crop == "rye_grass") {
    return(as.numeric(as.POSIXct(paste(year, "-11-05", sep = ""), format = "%Y-%m-%d")))
  } else if (crop == "rapeseed") {
    return(as.numeric(as.POSIXct(paste(year, "-09-15", sep = ""), format = "%Y-%m-%d")))
  } else if (crop == "potato") {
    return(as.numeric(as.POSIXct(paste(year, "-03-01", sep = ""), format = "%Y-%m-%d")))
  } else {
    return(NA)
  }
}

# Apply the function to create the crop_start_date column for each row in lhs_samples
# Convert lhs_samples to a data frame
lhs_df <- as.data.frame(lhs_samples)

# Name the columns
names(lhs_df) <- c("Value", "Year", "Crop")

# Calculate crop start dates based on lhs_df
crop_start_dates <- mapply(calculate_crop_start_date, lhs_df$Year, lhs_df$Crop)

#crop_start_dates <- mapply(calculate_crop_start_date, combined_df$year, combined_df$crop)


# Define a function to calculate crop start date and boundaries
# Define a function to calculate crop start date and boundaries
calculate_crop_info <- function(value, year, crop) {
  crop_start_date <- calculate_crop_start_date(year, crop)
  if (crop == "fallow") {
    min_boundary_crop <- as.numeric(crop_start_date) - 5 * 24 * 60 * 60  # 5 days earlier
    max_boundary_crop <- as.numeric(crop_start_date) + 5 * 24 * 60 * 60  # 5 days later
  } else {
    min_boundary_crop <- as.numeric(crop_start_date) - 60 * 24 * 60 * 60  # 60 days earlier
    max_boundary_crop <- as.numeric(crop_start_date) + 30 * 24 * 60 * 60  # 30 days later
  }
  return(list(crop_start_date = crop_start_date, min_boundary_crop = min_boundary_crop, max_boundary_crop = max_boundary_crop))
}

# Apply the function to each row of lhs_df
crop_info <- mapply(calculate_crop_info, lhs_df$Value, lhs_df$Year, lhs_df$Crop)


# Apply the function to each row of lhs_df
crop_info <- mapply(calculate_crop_info, lhs_df$Value, lhs_df$Year, lhs_df$Crop)

# Denormalize the Value column
# Create an empty list to store crop information
crop_info <- vector("list", length = nrow(lhs_df))

# Apply the function to each row of lhs_df and store the results in crop_info
for (i in seq_len(nrow(lhs_df))) {
  crop_info[[i]] <- calculate_crop_info(lhs_df[i, "Value"], lhs_df[i, "Year"], lhs_df[i, "Crop"])
}

# Denormalize the Value column
denormalized_values <- sapply(seq_len(nrow(lhs_df)), function(i) {
  value <- as.numeric(lhs_df[i, "Value"])
  min_boundary <- crop_info[[i]]$min_boundary_crop
  max_boundary <- crop_info[[i]]$max_boundary_crop
  return(value * (max_boundary - min_boundary) + min_boundary)
})

# Add the denormalized values to lhs_df
data_f$crop_start_date <- denormalized_values
timestamp1 <- as.POSIXct(data_f$crop_start_date, origin = "1970-01-01")
date_string1 <- format(timestamp1, format = "%d/%m/%Y")
cropstart=cbind(data_f$iteration,data_f$crop,date_string1)
#data_f$crop_start_date <- date_string1


# Initialize a vector to store the first appearance dates for each crop within each iteration
first_appearance <- vector("list", length = nrow(cropstart))

# Loop through each iteration
for (i in unique(cropstart[, 1])) {
  # Get the subset of cropstart for the current iteration
  iteration_cropstart <- cropstart[cropstart[, 1] == i, ]
  
  # Loop through each crop in the current iteration
  for (crop in unique(iteration_cropstart[, 2])) {
    # Find the first appearance of the current crop
    first_appearance_date <- iteration_cropstart[iteration_cropstart[, 2] == crop, 3][1]
    
    # Store the first appearance date in the vector
    first_appearance[[i]][[crop]] <- first_appearance_date
  }
}

# Replace the month and day in cropstart based on the first appearance of each crop within each iteration
for (i in 1:nrow(cropstart)) {
  iteration <- cropstart[i, 1]
  crop <- cropstart[i, 2]
  
  # Extract the year from the original date
  original_date <- as.Date(cropstart[i, 3], format = "%d/%m/%Y")
  year <- format(original_date, "%Y")
  
  # Use the month and day from the first appearance date, keeping the year fixed
  first_appearance_date <- as.Date(first_appearance[[iteration]][[crop]], format = "%d/%m/%Y")
  updated_date <- as.Date(paste(format(first_appearance_date, "%d/%m"), year, sep = "/"), format = "%d/%m/%Y")
  
  # Update the date in cropstart
  cropstart[i, 3] <- format(updated_date, "%d/%m/%Y")
}

# Display the updated cropstart dataframe
cropstart

data_f$crop_start_date <-  cropstart[,3]


set.seed(12345)
lhs_samples1 <- t(maximinLHS(1, nrow(data_f)))
lhs_samples1=cbind(lhs_samples1, data_f$year, data_f$crop)
calculate_dates_npk1 <- function(year, crop) {
  if (crop == "wheat") {
    year <- as.numeric(year) + 1
    npk_t1 <- as.numeric(as.POSIXct(paste(year, "-02-20", sep = ""), format = "%Y-%m-%d"))
  } else if (crop == "fallow") {
    npk_t1 <- 0
  } else if (crop == "maize") {
    npk_t1 <- as.numeric(as.POSIXct(paste(year, "-04-15", sep = ""), format = "%Y-%m-%d"))
  } else if (crop == "barley") {
    npk_t1 <- as.numeric(as.POSIXct(paste(year, "-03-20", sep = ""), format = "%Y-%m-%d"))
  } else if (crop == "rye_grass") {
    year <- as.numeric(year) + 1
    npk_t1 <- as.numeric(as.POSIXct(paste(year, "-02-20", sep = ""), format = "%Y-%m-%d"))
  } else if (crop == "rapeseed") {
    npk_t1 <- as.numeric(as.POSIXct(paste(year, "-12-25", sep = ""), format = "%Y-%m-%d"))
  } else if (crop == "potato") {
    npk_t1 <- as.numeric(as.POSIXct(paste(year, "-04-05", sep = ""), format = "%Y-%m-%d"))
  } else {
    npk_t1 <- NA
  }
  return(c(npk_t1))
}


# Apply the function to create the crop_start_date column for each row in lhs_samples
# Convert lhs_samples to a data frame
lhs_df1 <- as.data.frame(lhs_samples1)

# Name the columns
names(lhs_df1) <- c("Value", "Year", "Crop")

# Calculate crop start dates based on lhs_df
#crop_start_dates1 <- mapply(calculate_dates_npk1 , lhs_df$Year, lhs_df$Crop)

# Calculate crop start dates based on combined_df
crop_start_dates1 <- mapply(calculate_dates_npk1, data_f$year, data_f$crop)

# Define a function to calculate crop start date and boundaries
calculate_crop_info1 <- function(value, year, crop) {
  crop_start_date1 <- calculate_dates_npk1(year, crop)
  min_boundary_crop1 <- as.numeric(crop_start_date1) - 5 * 24 * 60 * 60
  max_boundary_crop1 <- as.numeric(crop_start_date1) + 30 * 24 * 60 * 60
  return(list(crop_start_date1 = crop_start_date1, min_boundary_crop1 = min_boundary_crop1, max_boundary_crop1 = max_boundary_crop1))
}

# Create an empty list to store crop information
crop_info <- vector("list", length = nrow(lhs_df1))

# Apply the function to each row of combined_df and store the results in crop_info
for (i in seq_len(nrow(lhs_df1))) {
  crop_info[[i]] <- calculate_crop_info1(lhs_df1[i, "Value"], lhs_df1[i, "Year"], lhs_df1[i, "Crop"])
}

# Denormalize the Value column
denormalized_values1 <- sapply(seq_len(nrow(lhs_df1)), function(i) {
  value <- as.numeric(lhs_df1[i, "Value"])
  min_boundary <- crop_info[[i]]$min_boundary_crop1
  max_boundary <- crop_info[[i]]$max_boundary_crop1
  if (lhs_df[i, "Crop"] == "fallow") {
    return(NA)
  } else {
    return(value * (max_boundary - min_boundary) + min_boundary)
  }
})

# Add the denormalized values to lhs_df1
npk_t1 <- denormalized_values1
##############################################################################
calculate_dates <- function(year, crop, denormalized_value) {
  if (crop == "wheat") {
    year <- year + 1
    npk_t1 <- denormalized_value
    npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
    npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
    npk_t4 <- NA
  } else if (crop == "fallow") {
    npk_t1 <- 0
    npk_t2 <- 0
    npk_t3 <- 0
    npk_t4 <- NA
  } else if (crop == "maize") {
    npk_t1 <- denormalized_value
    npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
    npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
    npk_t4 <- NA
  } else if (crop == "barley") {
    npk_t1 <- denormalized_value
    npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 30 days later
    npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 25 days later
    npk_t4 <- NA
  } else if (crop == "rye_grass") {
    year <- year + 1
    npk_t1 <- denormalized_value
    npk_t2 <- npk_t1 + (20 * 24 * 60 * 60)  # 20 days later
    npk_t3 <- npk_t2 + (20 * 24 * 60 * 60)  # 20 days later
    npk_t4 <- NA
  } else if (crop == "rapeseed") {
    npk_t1 <- denormalized_value
    npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 60 days later
    npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 30 days later
    npk_t4 <- npk_t3 + (10 * 24 * 60 * 60)  # 30 days later
  } else if (crop == "potato") {
    npk_t1 <- denormalized_value
    npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 20 days later
    npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 20 days later
    npk_t4 <- npk_t3 + (10 * 24 * 60 * 60)  # 20 days later
  } else {
    npk_t1 <- NA
    npk_t2 <- NA
    npk_t3 <- NA
    npk_t4 <- NA
  }
  
  return(c(npk_t1, npk_t2, npk_t3, npk_t4))
}
dates <- t(mapply(calculate_dates, data_f$year, data_f$crop, npk_t1))
colnames(dates) <- c("NPK_T1", "NPK_T2", "NPK_T3", "NPK_T4")
data_df <- data.frame(dates)

# Format and replace 01/01/1970 with empty string
for (i in 1:ncol(data_df)) {
  timestamp <- as.POSIXct(data_df[, i], origin = "1970-01-01")
  date_string <- format(timestamp, format = "%d/%m/%Y")
  data_df[, i] <- ifelse(date_string == "01/01/1970", "", date_string)
}
data_df[is.na(data_df)] <- ""
data_f <- cbind(data_f, data_df)  # Combine the rotation data with NPK dates
data_f[is.na(data_f)] <- ""
#######################
# Function to fix values for N, P, and K corresponding to different crops
fix_values_for_crops <- function(df) {
  # Iterate through each row of the dataframe
  for (i in 1:nrow(df)) {
    # Check the crop type and set values accordingly
    if (df[i, "crop"] == "wheat") {
      # Set values for wheat
      df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(105, 52.5, 52.5, 125, 62.5, 62.5, 125, 62.5, 62.5)
      df[i, c("N_4", "P_4", "K_4")] <- NA
    } else if (df[i, "crop"] == "barley") {
      # Set values for barley
      df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(105, 52.5, 52.5, 125, 62.5, 62.5, 125, 62.5, 62.5)
      df[i, c("N_4", "P_4", "K_4")] <- NA
    } else if (df[i, "crop"] == "maize") {
      # Set values for maize
      df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(20,10,10, 25, 12.5, 12.5, 35, 17.5, 17.5)
      df[i, c("N_4", "P_4", "K_4")] <- NA
    } else if (df[i, "crop"] == "potato") {
      # Set values for potato
      df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(125, 41.66667, 41.66667, 41.66667, 72.5, 24.16667,24.16667,24.16667,177.5,59.16667,59.16667,59.16667)
    } else if (df[i, "crop"] == "rapeseed") {
      # Set values for rapeseed
      df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(125,41.66667, 41.66667, 41.66667,60,20,20,20,175,58.33333,58.33333,58.33333)
    }
    else if (df[i, "crop"] == "fallow") {
      # Set values for rapeseed
      df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(NA)
    }
    # Add conditions for other crops if needed
  }
  return(df)
}

# Apply the function to fix values for different crops
data_f <- fix_values_for_crops(data_f)
########################
data_final_m <- replace_empty_cells(data_f)
data_final_m[is.na(data_final_m)] <- ""

#return(data_final_m)

# combined dataframe to CSV file
output_file_m_csv <- "C:/Users/mh1176/OneDrive - University of Exeter/Documents/wofost-rotations-multilonlat-itera.csv"
write.csv(data_final_m, output_file_m_csv, row.names = FALSE)
testfile_m <- read.csv(output_file_m_csv )
return(testfile_m)

}
lon_lat_pairs <- data.frame(
  lon = c(-1.83806), # Example longitudes
  lat = c(55.10674) # Example latitudes
)
rotations <- list(
  rotation_1 = c("wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat"))

# Call the function with rotations and lon_lat_pairs
generate_wofost_rotations_multi(rotations, lon_lat_pairs)




#############################################

########################################
# Define the number of lon-lat pairs
num_pairs <- 2

# Define the number of sets per pair
num_sets <- 5
it=num_pairs*num_sets
# Create a vector to store the iterations for all lon-lat pairs
iterations <- 1:it
iterations
# Create a dataframe with lon-lat pairs and iteration
lon_lat_pairs <- data.frame(
  lon = rep(c(-3.5, -3.52), each = num_sets),  # Replicate lon coordinates for each pair
  lat = rep(c(50.73, 50.73), each = num_sets),  # Replicate lat coordinates for each pair
  iteration = iterations # Sequence of iterations for all pairs
)

# Display the lon-lat pairs dataframe
print(lon_lat_pairs)

lon_lat_pairs <- data.frame(
  lon =c(-3.5),  # Replicate lon coordinates for each pair
  lat = c(50.73) # Sequence of iterations for all pairs
)


k=num_sets+1
p=2*num_sets
a=rep(k:p, each =  length(combined_df$iteration ))


# Define the number of lon-lat pairs
num_pairs <- 2

# Define the number of sets per pair
num_sets <- 5
ps=num_sets *num_sets
# Create an empty vector to store all iteration numbers
all_iterations <- c()

# Loop through each lon-lat pair
for (i in 1:num_pairs) {
  # Calculate the starting and ending iteration numbers for the current pair
  start_iteration <- (i - 1) * num_sets + 1
  end_iteration <- i * num_sets
  
  # Create a sequence of iteration numbers for the current pair
  pair_iterations <- seq(start_iteration, end_iteration)
  
  # Append the sequence to the vector of all iteration numbers
  all_iterations <- c(all_iterations, pair_iterations)
}

# Repeat the iteration numbers for each entry in combined_df$iteration
all_iterations <- rep(all_iterations, each = length(combined_df$iteration))

# Create a dataframe with lon-lat pairs and iteration numbers
lon_lat_pairs <- data.frame(
  lon = rep(c(-3.5, -3.52), each = ps),  # Replicate lon coordinates for each pair
  lat = rep(c(50.73, 50.71), each = ps),  # Replicate lat coordinates for each pair
  iteration = all_iterations # Sequence of iterations for all pairs
)

# Display the resulting dataframe
print(lon_lat_pairs)
data_f <- rbind(data_f, data_f)
data_f$lon=lon_lat_pairs$lon


# Number of times to replicate data_f
num_replications <- 1000

# Replicate data_f using replicate function
replicated_data <- replicate(num_replications, data_f, simplify = FALSE)

# Combine all replicated data frames into a single data frame
final_data_f <- do.call(rbind, replicated_data)

final_data_f 



# Define the number of lon-lat pairs
num_pairs <- 10

# Define the number of sets per pair
num_sets <- 5

# Create an empty vector to store all iteration numbers
all_iterations <- c()

# Loop through each lon-lat pair
for (i in 1:num_pairs) {
    # Calculate the starting and ending iteration numbers for the current pair
    start_iteration <- (i - 1) * num_sets + 1
    end_iteration <- i * num_sets
    
    # Create a sequence of iteration numbers for the current pair
    pair_iterations <- seq(start_iteration, end_iteration)
    
    # Append the sequence to the vector of all iteration numbers
    all_iterations <- c(all_iterations, pair_iterations)
}

# Repeat the iteration numbers for each entry in combined_df$iteration
all_iterations <- rep(all_iterations, each = length(combined_df$iteration))

# Display the resulting iteration vector
print(all_iterations)

####################################
generate_wofost_rotations_multi <- function(rotations,lon_lat_pairs) {
  # Function to adjust years for each rotation
  adjust_years <- function(rotation) {
    years <- rep(0, length(rotation))
    for (i in 1:length(rotation)) {
      if (i == 1) {
        years[i] <- 2020
      } else {
        if (rotation[i] == "fallow" && rotation[i - 1] %in% c("wheat", "rapeseed")) {
          years[i] <- years[i - 1] + 1  # Fallow in the next year after wheat or rapeseed
        } else if (rotation[i] == "fallow" && rotation[i - 1] %in% c("maize", "barley", "potato", "rye")) {
          years[i] <- years[i - 1]  # Fallow in the same year as starting year for maize, barley, potato, or rye
        } else if (rotation[i] == "fallow") {
          years[i] <- years[i - 1] + 1  # Fallow in the next year
        } else if (rotation[i] %in% c("barley", "maize", "potato", "rye")) {
          years[i] <- years[i - 1] + 1  # Crops after fallow in the same year
        } else {
          years[i] <- years[i - 1]  # Default: Same year
        }
      }
    }
    return(years)
  }
  
  # Apply the adjust_years function to each rotation
  rotation_years <- lapply(rotations, adjust_years)
  # Initialize a list to store rotation dataframes
  rotation_dfs <- list()
  
  # Iterate over rotations and lon-lat pairs
  for (rotation_name in names(rotations)) {
    for (i in 1:nrow(lon_lat_pairs)) {
      lon <- lon_lat_pairs$lon[i]
      lat <- lon_lat_pairs$lat[i]
      iteration <- i
      rotation <- rotations[[rotation_name]]
      years <- rotation_years[[rotation_name]]
      
      # Create a rotation dataframe for each lon-lat pair
      rotation_df <- data.frame(
        #rotation_id = paste0(rotation_name, "_", i, ".", seq_along(rotation)),
        rotation = rep(paste0(rotation_name, ".", i), length(rotation)),
        crop = rotation,
        year = years,
        lon = rep(lon, length(rotation)),
        lat = rep(lat, length(rotation)),
        iteration = rep(iteration, length(rotations[[rotation_name]]))
      )
      
      rotation_dfs[[paste0(rotation_name, "_", i)]] <- rotation_df
    }
  }
  
  combined_df <- do.call(rbind, rotation_dfs)
  combined_df <- combined_df[order(combined_df$iteration), ]
  
  # Define a function to add suffix to crop names
  add_variety <- function(crop) {
    if (crop == "wheat") {
      return("Winter_wheat_101")
    } else if (crop == "maize") {
      return("Grain_maize_201")
    } else if (crop == "barley") {
      return("Spring_barley_301")
    } else if (crop == "rye_grass") {
      return("Northern_RyeGrass")
    } else if (crop == "rapeseed") {
      return("Oilseed_rape_1001")
    } else if (crop == "potato") {
      return("Potato_701")
    } else if (crop == "fallow") {
      return("fallow")
    } else {
      return(NA)
    }
  }
  
  # Apply the function to create a new column with variety names
  combined_df$variety <- sapply(combined_df$crop, add_variety)
  # Define column names for the LHS inputs
  cols <- c("WAV", "SMLIM", "NAVAILI", "PAVAILI", "KAVAILI",
            "N_1", "N_2", "N_3", "N_4", "P_1", "P_2", "P_3","P_4", "K_1", "K_2", "K_3","K_4")
  
  # Define minimum and maximum values for each input variable
  min_values <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0)
  max_values <- c(60, 0.5, 100, 50, 100, 250, 250, 250, 250, 250, 250, 250, 250, 250,250,250,250)
  
  # Create a function to descale the inputs
  descale_inputs <- function(input_n, min_values, max_values) {
    if (!is.matrix(input_n)){
      input_n<- t(input_n)
    }
    num_rows <- ifelse(is.null(nrow(input_n)), 1, nrow(input_n)) 
    num_cols <- ifelse(is.null(ncol(input_n)), length(input_n), ncol(input_n)) 
    descaled_inputs <- matrix(0, nrow = num_rows, ncol = num_cols)
    
    for (i in 1:num_rows) {
      for (j in 1:num_cols) {
        descaled_inputs[i, j] <- (input_n[i, j] * (max_values[j] - min_values[j])) + min_values[j]
      }
    }
    
    return(descaled_inputs)
  }
  library(lhs)
  # Define the number of sets to generate
  num_sets <- 110
  
  # Create an empty list to store the LHS samples
  lhs_samplesa <- list()
  descaled_samples <- list()
  # Loop to generate LHS samples
  for (i in 1:num_sets) {
    #df3 <- maximinLHS(nrow(data_combined), length(cols))
    df3=maximinLHS(nrow(combined_df), length(cols))
    lhs_samplesa[[i]] <- df3
    descaled_samples[[i]] <- descale_inputs(lhs_samplesa[[i]], min_values, max_values)
  }
  
  #df3 <- maximinLHS(nrow(data_combined), 17) 
  input_data <- do.call(rbind, descaled_samples)
  colnames(input_data) <- cols
  
  # Function to replace empty cells for N_1 to K_4 for the crop "fallow"
  replace_empty_cells <- function(df) {
    # Iterate through each row of the dataframe
    for (i in 1:nrow(df)) {
      # Check if the crop is "fallow"
      if (df[i, "crop"] == "fallow") {
        # For "fallow", clear all N_1 to K_4
        for (j in 1:4) {
          df[i, paste0("N_", j)] <- ""
          df[i, paste0("P_", j)] <- ""
          df[i, paste0("K_", j)] <- ""
        }
      } else if (df[i, "crop"] %in% c("wheat", "barley", "maize")) {
        # For "wheat", "barley", and "maize", only clear N_4, P_4, and K_4
        df[i, "N_4"] <- ""
        df[i, "P_4"] <- ""
        df[i, "K_4"] <- ""
      }
    }
    return(df)
  }
  
  # Combine LHS inputs with the combined dataframe
  replicated_matrices=lapply(1:num_sets, function(x) combined_df)
  data_combined <- do.call(rbind, replicated_matrices)
  data_combined$iteration<- rep(1:num_sets, each =  length(combined_df$iteration ))
  data_f <- cbind(input_data, data_combined)
  
  
  
  # Generate LHS samples
  set.seed(123456)
  lhs_samples <- t(maximinLHS(1, nrow(data_f)))
  lhs_samples=cbind(lhs_samples, data_f$year, data_f$crop)
  
  
  calculate_crop_start_date <- function(year, crop) {
    if (crop == "wheat") {
      return(as.numeric(as.POSIXct(paste(year, "-11-05", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "fallow") {
      return(as.numeric(as.POSIXct(paste(year, "-09-01", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "maize") {
      return(as.numeric(as.POSIXct(paste(year, "-04-01", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "barley") {
      return(as.numeric(as.POSIXct(paste(year, "-02-15", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "rye_grass") {
      return(as.numeric(as.POSIXct(paste(year, "-11-05", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "rapeseed") {
      return(as.numeric(as.POSIXct(paste(year, "-09-15", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "potato") {
      return(as.numeric(as.POSIXct(paste(year, "-04-01", sep = ""), format = "%Y-%m-%d")))
    } else {
      return(NA)
    }
  }
  
  # Apply the function to create the crop_start_date column for each row in lhs_samples
  # Convert lhs_samples to a data frame
  lhs_df <- as.data.frame(lhs_samples)
  
  # Name the columns
  names(lhs_df) <- c("Value", "Year", "Crop")
  
  # Calculate crop start dates based on lhs_df
  crop_start_dates <- mapply(calculate_crop_start_date, lhs_df$Year, lhs_df$Crop)
  
  #crop_start_dates <- mapply(calculate_crop_start_date, combined_df$year, combined_df$crop)
  
  
  # Define a function to calculate crop start date and boundaries
  calculate_crop_info <- function(value, year, crop) {
    crop_start_date <- calculate_crop_start_date(year, crop)
    min_boundary_crop <- as.numeric(crop_start_date) - 5 * 24 * 60 * 60
    max_boundary_crop <- as.numeric(crop_start_date) + 5 * 24 * 60 * 60
    return(list(crop_start_date = crop_start_date, min_boundary_crop = min_boundary_crop, max_boundary_crop = max_boundary_crop))
  }
  
  
  
  # Apply the function to each row of lhs_df
  crop_info <- mapply(calculate_crop_info, lhs_df$Value, lhs_df$Year, lhs_df$Crop)
  
  # Denormalize the Value column
  # Create an empty list to store crop information
  crop_info <- vector("list", length = nrow(lhs_df))
  
  # Apply the function to each row of lhs_df and store the results in crop_info
  for (i in seq_len(nrow(lhs_df))) {
    crop_info[[i]] <- calculate_crop_info(lhs_df[i, "Value"], lhs_df[i, "Year"], lhs_df[i, "Crop"])
  }
  
  # Denormalize the Value column
  denormalized_values <- sapply(seq_len(nrow(lhs_df)), function(i) {
    value <- as.numeric(lhs_df[i, "Value"])
    min_boundary <- crop_info[[i]]$min_boundary_crop
    max_boundary <- crop_info[[i]]$max_boundary_crop
    return(value * (max_boundary - min_boundary) + min_boundary)
  })
  
  # Add the denormalized values to lhs_df
  data_f$crop_start_date <- denormalized_values
  timestamp1 <- as.POSIXct(data_f$crop_start_date, origin = "1970-01-01")
  date_string1 <- format(timestamp1, format = "%d/%m/%Y")
  cropstart=cbind(data_f$iteration,data_f$crop,date_string1)
  #data_f$crop_start_date <- date_string1
  
  
  # Initialize a vector to store the first appearance dates for each crop within each iteration
  first_appearance <- vector("list", length = nrow(cropstart))
  
  # Loop through each iteration
  for (i in unique(cropstart[, 1])) {
    # Get the subset of cropstart for the current iteration
    iteration_cropstart <- cropstart[cropstart[, 1] == i, ]
    
    # Loop through each crop in the current iteration
    for (crop in unique(iteration_cropstart[, 2])) {
      # Find the first appearance of the current crop
      first_appearance_date <- iteration_cropstart[iteration_cropstart[, 2] == crop, 3][1]
      
      # Store the first appearance date in the vector
      first_appearance[[i]][[crop]] <- first_appearance_date
    }
  }
  
  # Replace the month and day in cropstart based on the first appearance of each crop within each iteration
  for (i in 1:nrow(cropstart)) {
    iteration <- cropstart[i, 1]
    crop <- cropstart[i, 2]
    
    # Extract the year from the original date
    original_date <- as.Date(cropstart[i, 3], format = "%d/%m/%Y")
    year <- format(original_date, "%Y")
    
    # Use the month and day from the first appearance date, keeping the year fixed
    first_appearance_date <- as.Date(first_appearance[[iteration]][[crop]], format = "%d/%m/%Y")
    updated_date <- as.Date(paste(format(first_appearance_date, "%d/%m"), year, sep = "/"), format = "%d/%m/%Y")
    
    # Update the date in cropstart
    cropstart[i, 3] <- format(updated_date, "%d/%m/%Y")
  }
  
  # Display the updated cropstart dataframe
  cropstart
  
  data_f$crop_start_date <-  cropstart[,3]
  
  
  
  set.seed(12345)
  lhs_samples1 <- t(maximinLHS(1, nrow(data_f)))
  lhs_samples1=cbind(lhs_samples1, data_f$year, data_f$crop)
  calculate_dates_npk1 <- function(year, crop) {
    if (crop == "wheat") {
      year <- as.numeric(year) + 1
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-02-20", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "fallow") {
      npk_t1 <- 0
    } else if (crop == "maize") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-04-15", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "barley") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-03-20", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "rye_grass") {
      year <- as.numeric(year) + 1
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-02-20", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "rapeseed") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-12-05", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "potato") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-05-01", sep = ""), format = "%Y-%m-%d"))
    } else {
      npk_t1 <- NA
    }
    return(c(npk_t1))
  }
  
  
  # Apply the function to create the crop_start_date column for each row in lhs_samples
  # Convert lhs_samples to a data frame
  lhs_df1 <- as.data.frame(lhs_samples1)
  
  # Name the columns
  names(lhs_df1) <- c("Value", "Year", "Crop")
  
  # Calculate crop start dates based on lhs_df
  #crop_start_dates1 <- mapply(calculate_dates_npk1 , lhs_df$Year, lhs_df$Crop)
  
  # Calculate crop start dates based on combined_df
  crop_start_dates1 <- mapply(calculate_dates_npk1, data_f$year, data_f$crop)
  
  # Define a function to calculate crop start date and boundaries
  calculate_crop_info1 <- function(value, year, crop) {
    crop_start_date1 <- calculate_dates_npk1(year, crop)
    min_boundary_crop1 <- as.numeric(crop_start_date1) - 10 * 24 * 60 * 60
    max_boundary_crop1 <- as.numeric(crop_start_date1) + 30 * 24 * 60 * 60
    return(list(crop_start_date1 = crop_start_date1, min_boundary_crop1 = min_boundary_crop1, max_boundary_crop1 = max_boundary_crop1))
  }
  
  # Create an empty list to store crop information
  crop_info <- vector("list", length = nrow(lhs_df1))
  
  # Apply the function to each row of combined_df and store the results in crop_info
  for (i in seq_len(nrow(lhs_df1))) {
    crop_info[[i]] <- calculate_crop_info1(lhs_df1[i, "Value"], lhs_df1[i, "Year"], lhs_df1[i, "Crop"])
  }
  
  # Denormalize the Value column
  denormalized_values1 <- sapply(seq_len(nrow(lhs_df1)), function(i) {
    value <- as.numeric(lhs_df1[i, "Value"])
    min_boundary <- crop_info[[i]]$min_boundary_crop1
    max_boundary <- crop_info[[i]]$max_boundary_crop1
    if (lhs_df[i, "Crop"] == "fallow") {
      return(NA)
    } else {
      return(value * (max_boundary - min_boundary) + min_boundary)
    }
  })
  
  # Add the denormalized values to lhs_df1
  npk_t1 <- denormalized_values1
  ##############################################################################
  calculate_dates <- function(year, crop, denormalized_value) {
    if (crop == "wheat") {
      year <- year + 1
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- NA
    } else if (crop == "fallow") {
      npk_t1 <- 0
      npk_t2 <- 0
      npk_t3 <- 0
      npk_t4 <- NA
    } else if (crop == "maize") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- NA
    } else if (crop == "barley") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 30 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 25 days later
      npk_t4 <- NA
    } else if (crop == "rye_grass") {
      year <- year + 1
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 20 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 20 days later
      npk_t4 <- NA
    } else if (crop == "rapeseed") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- npk_t3 + (10 * 24 * 60 * 60)  # 10 days later
    } else if (crop == "potato") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- npk_t3 + (10 * 24 * 60 * 60)  # 10 days later
    } else {
      npk_t1 <- NA
      npk_t2 <- NA
      npk_t3 <- NA
      npk_t4 <- NA
    }
    
    return(c(npk_t1, npk_t2, npk_t3, npk_t4))
  }
  dates <- t(mapply(calculate_dates, data_f$year, data_f$crop, npk_t1))
  colnames(dates) <- c("NPK_T1", "NPK_T2", "NPK_T3", "NPK_T4")
  data_df <- data.frame(dates)
  
  # Format and replace 01/01/1970 with empty string
  for (i in 1:ncol(data_df)) {
    timestamp <- as.POSIXct(data_df[, i], origin = "1970-01-01")
    date_string <- format(timestamp, format = "%d/%m/%Y")
    data_df[, i] <- ifelse(date_string == "01/01/1970", "", date_string)
  }
  data_df[is.na(data_df)] <- ""
  data_f <- cbind(data_f, data_df)  # Combine the rotation data with NPK dates
  data_f[is.na(data_f)] <- ""
  #######################
  # Function to fix values for N, P, and K corresponding to different crops
  fix_values_for_crops <- function(df) {
    # Iterate through each row of the dataframe
    for (i in 1:nrow(df)) {
      # Check the crop type and set values accordingly
      if (df[i, "crop"] == "wheat") {
        # Set values for wheat
        df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(105, 52.5, 52.5, 125, 62.5, 62.5, 125, 62.5, 62.5)
        df[i, c("N_4", "P_4", "K_4")] <- NA
      } else if (df[i, "crop"] == "barley") {
        # Set values for barley
        df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(105, 52.5, 52.5, 125, 62.5, 62.5, 125, 62.5, 62.5)
        df[i, c("N_4", "P_4", "K_4")] <- NA
      } else if (df[i, "crop"] == "maize") {
        # Set values for maize
        df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(20,10,10, 25, 12.5, 12.5, 35, 17.5, 17.5)
        df[i, c("N_4", "P_4", "K_4")] <- NA
      } else if (df[i, "crop"] == "potato") {
        # Set values for potato
        df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(125, 41.66667, 41.66667, 41.66667, 72.5, 24.16667,24.16667,24.16667,177.5,59.16667,59.16667,59.16667)
      } else if (df[i, "crop"] == "rapeseed") {
        # Set values for rapeseed
        df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(125,41.66667, 41.66667, 41.66667,60,20,20,20,175,58.33333,58.33333,58.33333)
      }
      else if (df[i, "crop"] == "fallow") {
        # Set values for rapeseed
        df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(NA)
      }
      # Add conditions for other crops if needed
    }
    return(df)
  }
  
  # Apply the function to fix values for different crops
  data_f <- fix_values_for_crops(data_f)
  ########################
  data_final_m <- replace_empty_cells(data_f)
  data_final_m[is.na(data_final_m)] <- ""
  
  #return(data_final_m)
  
  # combined dataframe to CSV file
  output_file_m_csv <- "C:/Users/mh1176/OneDrive - University of Exeter/Documents/wofost-rotations-multilonlat-iter.csv"
  write.csv(data_final_m, output_file_m_csv, row.names = FALSE)
  testfile_m <- read.csv(output_file_m_csv )
  return(testfile_m)
  
}

rotations <- list(
  rotation_1 = c("wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat"))

# Call the function with rotations and lon_lat_pairs
generate_wofost_rotations_multi(rotations, lon_lat_pairs)


##########################wave2
generate_wofost_rotations_multi <- function(rotations,lon_lat_pairs) {
  # Function to adjust years for each rotation
  adjust_years <- function(rotation) {
    years <- rep(0, length(rotation))
    for (i in 1:length(rotation)) {
      if (i == 1) {
        years[i] <- 2020
      } else {
        if (rotation[i] == "fallow" && rotation[i - 1] %in% c("wheat", "rapeseed")) {
          years[i] <- years[i - 1] + 1  # Fallow in the next year after wheat or rapeseed
        } else if (rotation[i] == "fallow" && rotation[i - 1] %in% c("maize", "barley", "potato", "rye")) {
          years[i] <- years[i - 1]  # Fallow in the same year as starting year for maize, barley, potato, or rye
        } else if (rotation[i] == "fallow") {
          years[i] <- years[i - 1] + 1  # Fallow in the next year
        } else if (rotation[i] %in% c("barley", "maize", "potato", "rye")) {
          years[i] <- years[i - 1] + 1  # Crops after fallow in the same year
        } else {
          years[i] <- years[i - 1]  # Default: Same year
        }
      }
    }
    return(years)
  }
  
  # Apply the adjust_years function to each rotation
  rotation_years <- lapply(rotations, adjust_years)
  # Initialize a list to store rotation dataframes
  rotation_dfs <- list()
  
  # Iterate over rotations and lon-lat pairs
  for (rotation_name in names(rotations)) {
    for (i in 1:nrow(lon_lat_pairs)) {
      lon <- lon_lat_pairs$lon[i]
      lat <- lon_lat_pairs$lat[i]
      iteration <- i
      rotation <- rotations[[rotation_name]]
      years <- rotation_years[[rotation_name]]
      
      # Create a rotation dataframe for each lon-lat pair
      rotation_df <- data.frame(
        #rotation_id = paste0(rotation_name, "_", i, ".", seq_along(rotation)),
        rotation = rep(paste0(rotation_name, ".", i), length(rotation)),
        crop = rotation,
        year = years,
        lon = rep(lon, length(rotation)),
        lat = rep(lat, length(rotation)),
        iteration = rep(iteration, length(rotations[[rotation_name]]))
      )
      
      rotation_dfs[[paste0(rotation_name, "_", i)]] <- rotation_df
    }
  }
  
  combined_df <- do.call(rbind, rotation_dfs)
  combined_df <- combined_df[order(combined_df$iteration), ]
  
  # Define a function to add suffix to crop names
  add_variety <- function(crop) {
    if (crop == "wheat") {
      return("Winter_wheat_101")
    } else if (crop == "maize") {
      return("Grain_maize_201")
    } else if (crop == "barley") {
      return("Spring_barley_301")
    } else if (crop == "rye_grass") {
      return("Northern_RyeGrass")
    } else if (crop == "rapeseed") {
      return("Oilseed_rape_1001")
    } else if (crop == "potato") {
      return("Potato_701")
    } else if (crop == "fallow") {
      return("fallow")
    } else {
      return(NA)
    }
  }
  
  # Apply the function to create a new column with variety names
  combined_df$variety <- sapply(combined_df$crop, add_variety)
  # Define column names for the LHS inputs
  cols <- c("WAV", "SMLIM", "NAVAILI", "PAVAILI", "KAVAILI",
            "N_1", "N_2", "N_3", "N_4", "P_1", "P_2", "P_3","P_4", "K_1", "K_2", "K_3","K_4")
  
  # Define minimum and maximum values for each input variable
  min_values <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0)
  max_values <- c(60, 0.5, 100, 50, 100, 250, 250, 250, 250, 250, 250, 250, 250, 250,250,250,250)
  
  # Create a function to descale the inputs
  descale_inputs <- function(input_n, min_values, max_values) {
    if (!is.matrix(input_n)){
      input_n<- t(input_n)
    }
    num_rows <- ifelse(is.null(nrow(input_n)), 1, nrow(input_n)) 
    num_cols <- ifelse(is.null(ncol(input_n)), length(input_n), ncol(input_n)) 
    descaled_inputs <- matrix(0, nrow = num_rows, ncol = num_cols)
    
    for (i in 1:num_rows) {
      for (j in 1:num_cols) {
        descaled_inputs[i, j] <- (input_n[i, j] * (max_values[j] - min_values[j])) + min_values[j]
      }
    }
    
    return(descaled_inputs)
  }
  library(lhs)
  # Define the number of sets to generate
  num_sets <- 3
  
  # Create an empty list to store the LHS samples
  lhs_samplesa <- list()
  descaled_samples <- list()
  # Loop to generate LHS samples
  for (i in 1:num_sets) {
    df33=maximinLHS(nrow(combined_df), 12)
    df34=wave3_design[1:nrow(combined_df),1:5]
    df3=cbind(df34,df33)
    lhs_samplesa[[i]] <- df3
    descaled_samples[[i]] <- descale_inputs(lhs_samplesa[[i]], min_values, max_values)
  }
  
  #df3 <- maximinLHS(nrow(data_combined), 17) 
  input_data <- do.call(rbind, descaled_samples)
  colnames(input_data) <- cols
  
  # Function to replace empty cells for N_1 to K_4 for the crop "fallow"
  replace_empty_cells <- function(df) {
    # Iterate through each row of the dataframe
    for (i in 1:nrow(df)) {
      # Check if the crop is "fallow"
      if (df[i, "crop"] == "fallow") {
        # For "fallow", clear all N_1 to K_4
        for (j in 1:4) {
          df[i, paste0("N_", j)] <- ""
          df[i, paste0("P_", j)] <- ""
          df[i, paste0("K_", j)] <- ""
        }
      } else if (df[i, "crop"] %in% c("wheat", "barley", "maize")) {
        # For "wheat", "barley", and "maize", only clear N_4, P_4, and K_4
        df[i, "N_4"] <- ""
        df[i, "P_4"] <- ""
        df[i, "K_4"] <- ""
      }
    }
    return(df)
  }
  
  # Combine LHS inputs with the combined dataframe
  replicated_matrices=lapply(1:num_sets, function(x) combined_df)
  data_combined <- do.call(rbind, replicated_matrices)
  data_combined$iteration<- rep(1:num_sets, each =  length(combined_df$iteration ))
  data_f <- cbind(input_data, data_combined)
  
  
  
  # Generate LHS samples
  set.seed(123456)
  lhs_samples <- wave3_design[1:nrow(data_f),6]
  lhs_samples=cbind(lhs_samples, data_f$year, data_f$crop)
  
  
  calculate_crop_start_date <- function(year, crop) {
    if (crop == "wheat") {
      return(as.numeric(as.POSIXct(paste(year, "-11-05", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "fallow") {
      return(as.numeric(as.POSIXct(paste(year, "-09-01", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "maize") {
      return(as.numeric(as.POSIXct(paste(year, "-04-01", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "barley") {
      return(as.numeric(as.POSIXct(paste(year, "-02-15", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "rye_grass") {
      return(as.numeric(as.POSIXct(paste(year, "-11-05", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "rapeseed") {
      return(as.numeric(as.POSIXct(paste(year, "-09-15", sep = ""), format = "%Y-%m-%d")))
    } else if (crop == "potato") {
      return(as.numeric(as.POSIXct(paste(year, "-04-01", sep = ""), format = "%Y-%m-%d")))
    } else {
      return(NA)
    }
  }
  
  # Apply the function to create the crop_start_date column for each row in lhs_samples
  # Convert lhs_samples to a data frame
  lhs_df <- as.data.frame(lhs_samples)
  
  # Name the columns
  names(lhs_df) <- c("Value", "Year", "Crop")
  
  # Calculate crop start dates based on lhs_df
  crop_start_dates <- mapply(calculate_crop_start_date, lhs_df$Year, lhs_df$Crop)
  
  #crop_start_dates <- mapply(calculate_crop_start_date, combined_df$year, combined_df$crop)
  
  
  # Define a function to calculate crop start date and boundaries
  calculate_crop_info <- function(value, year, crop) {
    crop_start_date <- calculate_crop_start_date(year, crop)
    min_boundary_crop <- as.numeric(crop_start_date) - 5 * 24 * 60 * 60
    max_boundary_crop <- as.numeric(crop_start_date) + 5 * 24 * 60 * 60
    return(list(crop_start_date = crop_start_date, min_boundary_crop = min_boundary_crop, max_boundary_crop = max_boundary_crop))
  }
  
  
  
  # Apply the function to each row of lhs_df
  crop_info <- mapply(calculate_crop_info, lhs_df$Value, lhs_df$Year, lhs_df$Crop)
  
  # Denormalize the Value column
  # Create an empty list to store crop information
  crop_info <- vector("list", length = nrow(lhs_df))
  
  # Apply the function to each row of lhs_df and store the results in crop_info
  for (i in seq_len(nrow(lhs_df))) {
    crop_info[[i]] <- calculate_crop_info(lhs_df[i, "Value"], lhs_df[i, "Year"], lhs_df[i, "Crop"])
  }
  
  # Denormalize the Value column
  denormalized_values <- sapply(seq_len(nrow(lhs_df)), function(i) {
    value <- as.numeric(lhs_df[i, "Value"])
    min_boundary <- crop_info[[i]]$min_boundary_crop
    max_boundary <- crop_info[[i]]$max_boundary_crop
    return(value * (max_boundary - min_boundary) + min_boundary)
  })
  
  # Add the denormalized values to lhs_df
  data_f$crop_start_date <- denormalized_values
  timestamp1 <- as.POSIXct(data_f$crop_start_date, origin = "1970-01-01")
  date_string1 <- format(timestamp1, format = "%d/%m/%Y")
  cropstart=cbind(data_f$iteration,data_f$crop,date_string1)
  #data_f$crop_start_date <- date_string1
  
  
  # Initialize a vector to store the first appearance dates for each crop within each iteration
  first_appearance <- vector("list", length = nrow(cropstart))
  
  # Loop through each iteration
  for (i in unique(cropstart[, 1])) {
    # Get the subset of cropstart for the current iteration
    iteration_cropstart <- cropstart[cropstart[, 1] == i, ]
    
    # Loop through each crop in the current iteration
    for (crop in unique(iteration_cropstart[, 2])) {
      # Find the first appearance of the current crop
      first_appearance_date <- iteration_cropstart[iteration_cropstart[, 2] == crop, 3][1]
      
      # Store the first appearance date in the vector
      first_appearance[[i]][[crop]] <- first_appearance_date
    }
  }
  
  # Replace the month and day in cropstart based on the first appearance of each crop within each iteration
  for (i in 1:nrow(cropstart)) {
    iteration <- cropstart[i, 1]
    crop <- cropstart[i, 2]
    
    # Extract the year from the original date
    original_date <- as.Date(cropstart[i, 3], format = "%d/%m/%Y")
    year <- format(original_date, "%Y")
    
    # Use the month and day from the first appearance date, keeping the year fixed
    first_appearance_date <- as.Date(first_appearance[[iteration]][[crop]], format = "%d/%m/%Y")
    updated_date <- as.Date(paste(format(first_appearance_date, "%d/%m"), year, sep = "/"), format = "%d/%m/%Y")
    
    # Update the date in cropstart
    cropstart[i, 3] <- format(updated_date, "%d/%m/%Y")
  }
  
  # Display the updated cropstart dataframe
  cropstart
  
  data_f$crop_start_date <-  cropstart[,3]
  
  
  
  set.seed(12345)
  lhs_samples1 <- wave3_design[1:nrow(data_f),7]
  lhs_samples1=cbind(lhs_samples1, data_f$year, data_f$crop)
  calculate_dates_npk1 <- function(year, crop) {
    if (crop == "wheat") {
      year <- as.numeric(year) + 1
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-02-20", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "fallow") {
      npk_t1 <- 0
    } else if (crop == "maize") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-04-15", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "barley") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-03-20", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "rye_grass") {
      year <- as.numeric(year) + 1
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-02-20", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "rapeseed") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-12-05", sep = ""), format = "%Y-%m-%d"))
    } else if (crop == "potato") {
      npk_t1 <- as.numeric(as.POSIXct(paste(year, "-05-01", sep = ""), format = "%Y-%m-%d"))
    } else {
      npk_t1 <- NA
    }
    return(c(npk_t1))
  }
  
  
  # Apply the function to create the crop_start_date column for each row in lhs_samples
  # Convert lhs_samples to a data frame
  lhs_df1 <- as.data.frame(lhs_samples1)
  
  # Name the columns
  names(lhs_df1) <- c("Value", "Year", "Crop")
  
  # Calculate crop start dates based on lhs_df
  #crop_start_dates1 <- mapply(calculate_dates_npk1 , lhs_df$Year, lhs_df$Crop)
  
  # Calculate crop start dates based on combined_df
  crop_start_dates1 <- mapply(calculate_dates_npk1, data_f$year, data_f$crop)
  
  # Define a function to calculate crop start date and boundaries
  calculate_crop_info1 <- function(value, year, crop) {
    crop_start_date1 <- calculate_dates_npk1(year, crop)
    min_boundary_crop1 <- as.numeric(crop_start_date1) - 10 * 24 * 60 * 60
    max_boundary_crop1 <- as.numeric(crop_start_date1) + 30 * 24 * 60 * 60
    return(list(crop_start_date1 = crop_start_date1, min_boundary_crop1 = min_boundary_crop1, max_boundary_crop1 = max_boundary_crop1))
  }
  
  # Create an empty list to store crop information
  crop_info <- vector("list", length = nrow(lhs_df1))
  
  # Apply the function to each row of combined_df and store the results in crop_info
  for (i in seq_len(nrow(lhs_df1))) {
    crop_info[[i]] <- calculate_crop_info1(lhs_df1[i, "Value"], lhs_df1[i, "Year"], lhs_df1[i, "Crop"])
  }
  
  # Denormalize the Value column
  denormalized_values1 <- sapply(seq_len(nrow(lhs_df1)), function(i) {
    value <- as.numeric(lhs_df1[i, "Value"])
    min_boundary <- crop_info[[i]]$min_boundary_crop1
    max_boundary <- crop_info[[i]]$max_boundary_crop1
    if (lhs_df[i, "Crop"] == "fallow") {
      return(NA)
    } else {
      return(value * (max_boundary - min_boundary) + min_boundary)
    }
  })
  
  # Add the denormalized values to lhs_df1
  npk_t1 <- denormalized_values1
  ##############################################################################
  calculate_dates <- function(year, crop, denormalized_value) {
    if (crop == "wheat") {
      year <- year + 1
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- NA
    } else if (crop == "fallow") {
      npk_t1 <- 0
      npk_t2 <- 0
      npk_t3 <- 0
      npk_t4 <- NA
    } else if (crop == "maize") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- NA
    } else if (crop == "barley") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 30 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 25 days later
      npk_t4 <- NA
    } else if (crop == "rye_grass") {
      year <- year + 1
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 20 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 20 days later
      npk_t4 <- NA
    } else if (crop == "rapeseed") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- npk_t3 + (10 * 24 * 60 * 60)  # 10 days later
    } else if (crop == "potato") {
      npk_t1 <- denormalized_value
      npk_t2 <- npk_t1 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t3 <- npk_t2 + (10 * 24 * 60 * 60)  # 10 days later
      npk_t4 <- npk_t3 + (10 * 24 * 60 * 60)  # 10 days later
    } else {
      npk_t1 <- NA
      npk_t2 <- NA
      npk_t3 <- NA
      npk_t4 <- NA
    }
    
    return(c(npk_t1, npk_t2, npk_t3, npk_t4))
  }
  dates <- t(mapply(calculate_dates, data_f$year, data_f$crop, npk_t1))
  colnames(dates) <- c("NPK_T1", "NPK_T2", "NPK_T3", "NPK_T4")
  data_df <- data.frame(dates)
  
  # Format and replace 01/01/1970 with empty string
  for (i in 1:ncol(data_df)) {
    timestamp <- as.POSIXct(data_df[, i], origin = "1970-01-01")
    date_string <- format(timestamp, format = "%d/%m/%Y")
    data_df[, i] <- ifelse(date_string == "01/01/1970", "", date_string)
  }
  data_df[is.na(data_df)] <- ""
  data_f <- cbind(data_f, data_df)  # Combine the rotation data with NPK dates
  data_f[is.na(data_f)] <- ""
  #######################
  # Function to fix values for N, P, and K corresponding to different crops
  fix_values_for_crops <- function(df) {
    # Iterate through each row of the dataframe
    for (i in 1:nrow(df)) {
      # Check the crop type and set values accordingly
      if (df[i, "crop"] == "wheat") {
        # Set values for wheat
        df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(105, 52.5, 52.5, 125, 62.5, 62.5, 125, 62.5, 62.5)
        df[i, c("N_4", "P_4", "K_4")] <- NA
      } else if (df[i, "crop"] == "barley") {
        # Set values for barley
        df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(105, 52.5, 52.5, 125, 62.5, 62.5, 125, 62.5, 62.5)
        df[i, c("N_4", "P_4", "K_4")] <- NA
      } else if (df[i, "crop"] == "maize") {
        # Set values for maize
        df[i, c("N_1", "N_2", "N_3", "P_1", "P_2", "P_3", "K_1", "K_2", "K_3")] <- c(20,10,10, 25, 12.5, 12.5, 35, 17.5, 17.5)
        df[i, c("N_4", "P_4", "K_4")] <- NA
      } else if (df[i, "crop"] == "potato") {
        # Set values for potato
        df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(125, 41.66667, 41.66667, 41.66667, 72.5, 24.16667,24.16667,24.16667,177.5,59.16667,59.16667,59.16667)
      } else if (df[i, "crop"] == "rapeseed") {
        # Set values for rapeseed
        df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(125,41.66667, 41.66667, 41.66667,60,20,20,20,175,58.33333,58.33333,58.33333)
      }
      else if (df[i, "crop"] == "fallow") {
        # Set values for rapeseed
        df[i, c("N_1", "N_2", "N_3","N_4", "P_1", "P_2", "P_3", "P_4","K_1", "K_2", "K_3", "K_4")] <- c(NA)
      }
      # Add conditions for other crops if needed
    }
    return(df)
  }
  
  # Apply the function to fix values for different crops
  data_f <- fix_values_for_crops(data_f)
  ########################
  data_final_m <- replace_empty_cells(data_f)
  data_final_m[is.na(data_final_m)] <- ""
  
  #return(data_final_m)
  
  # combined dataframe to CSV file
  output_file_m_csv <- "C:/Users/mh1176/OneDrive - University of Exeter/Documents/wofost-rotations-multilonlat-iter.csv"
  write.csv(data_final_m, output_file_m_csv, row.names = FALSE)
  testfile_m <- read.csv(output_file_m_csv )
  return(testfile_m)
  
}

rotations <- list(
  rotation_1 = c("wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat",
                 "fallow","wheat", "fallow", "potato","fallow","wheat", "fallow", "rapeseed", "fallow", "wheat", "fallow","wheat"))

# Call the function with rotations and lon_lat_pairs
generate_wofost_rotations_multi(rotations, lon_lat_pairs)


