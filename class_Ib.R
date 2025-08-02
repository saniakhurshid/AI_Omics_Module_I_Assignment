getwd()

# create sub folders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# import patient data
patient_data = read.csv(file.choose())
str(patient_data)

# required changes
patient_data$gender = as.factor(patient_data$gender)
patient_data$diagnosis = as.factor(patient_data$diagnosis)
patient_data$smoker_binary = as.factor(ifelse(patient_data$smoker== "Yes",1,0))
str(patient_data)

# save clean data file
write.csv(patient_data, "clean_data/patient_info_clean.csv", row.names = FALSE)
