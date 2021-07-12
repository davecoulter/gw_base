from database_methods import *
upload_csv = "local_data/PanSTARRS_allData_v2_pjquinonez.csv"
upload_sql = """LOAD DATA LOCAL INFILE '%s' 
                    INTO TABLE PS1_Galaxy_v3 
                    FIELDS TERMINATED BY ',' 
                    LINES TERMINATED BY '\n' 
                    IGNORE 1 ROWS;"""

success = bulk_upload(upload_sql % upload_csv)
#(objID, uniquePspsOBid, raStack, decStack, raMean, decMean, ra, dec, ng, gMeanPSFMag, gMeanPSFMagErr, gMeanKronMag, gMeanKronMagErr, gMeanApMag, gMeanApMagErr, nr, rMeanPSFMag, rMeanPSFMagErr, rMeanKronMag, rMeanKronMagErr, rMeanApMag, rMeanApMagErr, ni, iMeanPSFMag, iMeanPSFMagErr, iMeanKronMag, iMeanKronMagErr, iMeanApMag, iMeanApMagErr, nz, zMeanPSFMag, zMeanPSFMagErr, zMeanKronMag, zMeanKronMagErr, zMeanApMag, zMeanApMagErr, ny, yMeanPSFMag, yMeanPSFMagErr, yMeanKronMag, yMeanKronMagErr, yMeanApMag, yMeanApMagErr, gQfPerfect, rQfPerfect, iQfPerfect, zQfPerfect, yQfPerfect, qualityFlag, objInfoFlag, gpetRadius, rpetRadius, ipetRadius, zpetRadius, ypetRadius, gpetR50, rpetR50, ipetR50, zpetR50, ypetR50, primaryDetection, bestDetection, gKronFlux, gKronFluxErr, rKronFlux, rKronFluxErr, iKronFlux, iKronFluxErr, zKronFlux, zKronFluxErr, yKronFlux, yKronFluxErr, class, prob_Galaxy, prob_Star, prob_QSO, z_phot, z_photErr, z_phot0, extrapolation_Photoz);