#### practice Metadata ####

# Load file
Metadata = read.csv(file.choose())

# copy file and structure
mydata <- Metadata
str(mydata)

# factor cols in a vector
factor_cols = c('height', 'gender')

# for Loop
for (col in factor_cols) {mydata[[col]] <- as.factor(mydata[[col]])}
str(mydata)

# binary cols in a vector
binary_cols = c('gender')

# ifelse in for Loop
for (col in binary_cols) {mydata[[col]] = ifelse(mydata[[col]] == 'Male',1,0)}
str(mydata)

# comparison
str(Metadata)
str(mydata)


#### practice patient_info data ####

# Load file
patient_info = read.csv(file.choose())

# copy file and structure
patient_data = patient_info
str(patient_data)

# factors cols in a vector
factor_cols2 = c('gender','diagnosis','smoker')

# for Loop
for (col in factor_cols2 ) {patient_data[[col]] <- as.factor(patient_data[[col]])}
str(patient_data)

# binary cols in a vector
binary_cols2 = c('gender','diagnosis','smoker')

# ifelse in for Loop
for (col in binary_cols2) {
  if (col == "gender") {
    patient_data[[col]] <- ifelse(patient_data[[col]] == "Male", 1, 0)
  } else if (col == "diagnosis") {
    patient_data[[col]] <- ifelse(patient_data[[col]] == "Cancer", 1, 0)
  } else if (col == "smoker") {
    patient_data[[col]] <- ifelse(patient_data[[col]] == "Yes", 1, 0)
  }
}

# Comparison
str(patient_info)
str(patient_data)



