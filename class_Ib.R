getwd()

# create sub folders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# import patient data
patient_data = read.csv(file.choose())
View(patient_data)
str(patient_data)

# required changes in variables
patient_data$gender = as.factor(patient_data$gender)
patient_data$diagnosis <- as.factor(patient_data$diagnosis)
patient_data$smoker_binary = as.factor(ifelse(patient_data$smoker == "yes", 1,0))
str(patient_data)

# clean data folder
write.csv(patient_data, "clean_data/patient_info.csv", row.names = FALSE)


