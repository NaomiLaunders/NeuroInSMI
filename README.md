# The temporal relationship between severe mental illness and neurological conditions in a UK primary care cohort

These files include the code lists and scripts used in the above study of CPRD data.

## Understanding the files
Scripts are organised into chunks: Data processing (scripts beginning 1), Creating the analysis data set (scripts beginning 2), and analysis (scripts beginning 3).

Code lists were generated in CPRD. They are specific to CPRD and to the exact database builds. When re-using codelists it's advisable to search for new terms in the system being used.

CPRD code lists use medcodes. These long numeric identifiers must be kept as character vectors in R, Stata and Excel to avoid truncation.

This project uses pseudonymised patient records and so no data is available for uploading.
