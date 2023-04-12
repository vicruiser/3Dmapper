############################################
#
#	Run test 
#
#

# Download PDBs with ID
#
# Create human proteome DB with blast
#
# Create protein structure db
#makeinterfacesdb 

##############################
# Split protein structure db #
##############################

makepsdb

##############################
# Split variants db #
##############################

makevariantsdb

##############################
# map
##############################
mapper 

