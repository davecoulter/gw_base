from database_methods import *


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--healpix_map_id', default="-1", type="int",
                          help='''The integer primary key for the map to remove.''')

        return (parser)

    def main(self):
        is_error = False

        if self.options.healpix_map_id < 0:
            is_error = True
            print("HealpixMap_id is required.")

        healpix_map_select = '''
            SELECT GWID, URL, Filename FROM HealpixMap WHERE id = '%s'
        '''

        healpix_pixel_select = '''
            SELECT
                id
            FROM
                HealpixPixel
            WHERE
                HealpixMap_id = %s
        '''

        healpix_completeness_select = '''
            SELECT 
                id 
            FROM 
                HealpixPixel_Completeness 
            WHERE HealpixPixel_id IN (%s)
        '''

        healpix_pixel_galaxy_distance_select = '''
            SELECT 
                id 
            FROM 
                HealpixPixel_GalaxyDistance2 
            WHERE 
                HealpixPixel_id IN (%s)
        '''

        healpix_pixel_galaxy_distance_weight_select = '''
            SELECT 
                id 
            FROM 
                HealpixPixel_GalaxyDistance2_Weight  
            WHERE 
                HealpixPixel_GalaxyDistance2_id IN (%s)
        '''

        observed_tile_healpix_pixel_select = '''
            SELECT 
                id 
            FROM 
                ObservedTile_HealpixPixel 
            WHERE 
                HealpixPixel_id IN (%s) 
        '''

        observed_tile_select = '''
            SELECT 
                id 
            FROM 
                ObservedTile  
            WHERE 
                HealpixMap_id = %s   
        '''

        static_tile_healpix_pixel_select = '''
            SELECT 
                id 
            FROM 
                StaticTile_HealpixPixel 
            WHERE 
                HealpixPixel_id IN (%s) 
        '''

        healpix_completeness_delete = '''
            DELETE FROM HealpixPixel_Completeness 
            WHERE id IN (%s) AND id > 0 
        '''

        healpix_pixel_galaxy_distance_weight_delete = '''
            DELETE FROM HealpixPixel_GalaxyDistance2_Weight 
            WHERE id IN (%s) AND id > 0 
        '''

        healpix_pixel_galaxy_distance_delete = '''
            DELETE FROM HealpixPixel_GalaxyDistance2 
            WHERE id IN (%s) AND id > 0 
        '''

        observed_tile_healpix_pixel_delete = '''
            DELETE FROM ObservedTile_HealpixPixel 
            WHERE id IN (%s) AND id > 0 
        '''

        observed_tile_delete = '''
            DELETE FROM ObservedTile  
            WHERE id IN (%s) AND id > 0 
        '''

        static_tile_healpix_pixel_delete = '''
            DELETE FROM StaticTile_HealpixPixel 
            WHERE id IN (%s) AND id > 0 
        '''

        healpix_pixel_delete = '''
            DELETE FROM HealpixPixel 
            WHERE id IN (%s) AND id > 0 
        '''

        healpix_map_delete = '''
            DELETE FROM HealpixMap    
            WHERE id = %s 
        '''

        if is_error:
            print("Exiting...")
            return 1


        observed_tile_result = []
        healpix_pixel_result = []
        healpix_pixel_id_string = ""
        healpix_map_result = query_db([healpix_map_select % self.options.healpix_map_id])[0]
        if len(healpix_map_result) > 0:
            healpix_pixel_result = query_db([healpix_pixel_select % self.options.healpix_map_id])[0]
            healpix_pixel_id_string = ",".join([str(hp[0]) for hp in healpix_pixel_result])
            observed_tile_result = query_db([observed_tile_select % self.options.healpix_map_id])[0]

        observed_tile_healpix_pixel_result = []
        static_tile_healpix_pixel_result = []
        healpix_completeness_result = []
        healpix_pixel_galaxy_distance_result = []
        if len(healpix_pixel_result) > 0:
            observed_tile_healpix_pixel_result = query_db([observed_tile_healpix_pixel_select % healpix_pixel_id_string])[0]
            static_tile_healpix_pixel_result = query_db([static_tile_healpix_pixel_select % healpix_pixel_id_string])[0]
            healpix_completeness_result = query_db([healpix_completeness_select % healpix_pixel_id_string])[0]
            healpix_pixel_galaxy_distance_result = query_db([healpix_pixel_galaxy_distance_select % healpix_pixel_id_string])[0]

        healpix_pixel_galaxy_distance_weight_result = []
        if len(healpix_pixel_galaxy_distance_result) > 0:
            healpix_pixel_galaxy_distance_id_string = ",".join([str(hpgd[0]) for hpgd in healpix_pixel_galaxy_distance_result])
            healpix_pixel_galaxy_distance_weight_result = query_db([healpix_pixel_galaxy_distance_weight_select % healpix_pixel_galaxy_distance_id_string])[0]


        if len(healpix_map_result) > 0:
            print("Are you sure you want to delete:\n\tGWID: %s\n\tURL: %s\n\tFile: %s" % healpix_map_result[0])

        print('''Number of other records to be deleted:
            HealpixPixel_Completeness: %s
            HealpixPixel_GalaxyDistance2_Weight: %s
            HealpixPixel_GalaxyDistance2: %s
            ObservedTile_HealpixPixel: %s
            ObservedTile: %s
            StaticTile_HealpixPixel: %s
            HealpixPixel: %s
                ''' % (len(healpix_completeness_result),
                       len(healpix_pixel_galaxy_distance_weight_result),
                       len(healpix_pixel_galaxy_distance_result),
                       len(observed_tile_healpix_pixel_result),
                       len(observed_tile_result),
                       len(static_tile_healpix_pixel_result),
                       len(healpix_pixel_result)))

        if len(healpix_completeness_result) > 0:
            print("DELETE %s HealpixPixel_Completeness" % len(healpix_completeness_result))
            delete_rows(healpix_completeness_delete, healpix_completeness_result)

        if len(healpix_pixel_galaxy_distance_weight_result) > 0:
            print("DELETE %s HealpixPixel_GalaxyDistance2_Weight" % len(healpix_pixel_galaxy_distance_weight_result))
            delete_rows(healpix_pixel_galaxy_distance_weight_delete, healpix_pixel_galaxy_distance_weight_result)

        if len(healpix_pixel_galaxy_distance_result) > 0:
            print("DELETE %s HealpixPixel_GalaxyDistance2" % len(healpix_pixel_galaxy_distance_result))
            delete_rows(healpix_pixel_galaxy_distance_delete, healpix_pixel_galaxy_distance_result)

        if len(observed_tile_healpix_pixel_result) > 0:
            print("DELETE %s ObservedTile_HealpixPixel" % len(observed_tile_healpix_pixel_result))
            delete_rows(observed_tile_healpix_pixel_delete, observed_tile_healpix_pixel_result)

        if len(observed_tile_result) > 0:
            print("DELETE %s ObservedTile" % len(observed_tile_result))
            delete_rows(observed_tile_delete, observed_tile_result)

        if len(static_tile_healpix_pixel_result) > 0:
            print("DELETE %s StaticTile_HealpixPixel" % len(static_tile_healpix_pixel_result))
            delete_rows(static_tile_healpix_pixel_delete, static_tile_healpix_pixel_result)

        if len(healpix_pixel_result) > 0:
            print("DELETE %s HealpixPixel" % len(healpix_pixel_result))
            delete_rows(healpix_pixel_delete, healpix_pixel_result)

        if len(healpix_map_result) > 0:
            print("DELETE HealpixMap id=%s" % str(self.options.healpix_map_id))
            delete_map = healpix_map_delete % str(self.options.healpix_map_id)
            query_db([delete_map], commit=True)


if __name__ == "__main__":
    useagestring = """python DeleteMap.py [options]
    
python DeleteMap.py --healpix_map_id <database id>
"""

    start = time.time()

    teglon = Teglon()
    parser = teglon.add_options(usage=useagestring)
    options, args = parser.parse_args()
    teglon.options = options

    teglon.main()

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `DeleteMap` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")


