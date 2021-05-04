from Pixel_Element import *
from SQL_Polygon import *
from shapely import geometry


class Detector(Teglon_Shape):

    # returns a collection of N=len(resolution) vertex polygon, in degrees, centered at the origin, and
    #   defined in clockwise order from the first-quadrant:
    #    [
    #      (ra_1, dec_1), (ra_2, dec_2), ... (ra_N, dec_N)
    #    ]
    @staticmethod
    def get_detector_vertices_circular(deg_radius, resolution=50):
        c1 = Point(0, 0).buffer(deg_radius)
        # clip the last 2 coords off to protect against degenerate corners
        return np.asarray([np.asarray([(c[0], c[1]) for c in c1.exterior.coords[:-2]])])

    # returns a collection of N=4 vertex polygon, in degreees, centered at the origin, and
    #   defined in clockwise order from the first-quadrant:
    #    [
    #      (ra_1, dec_1), (ra_2, dec_2), ... (ra_N, dec_N)
    #    ]
    @staticmethod
    def get_detector_vertices_rectangular(deg_width, deg_height):

        x_step = deg_width / 2.0
        y_step = deg_height / 2.0

        clockwise_square = [(x_step, y_step),
                            (x_step, -1.0 * y_step),
                            (-1.0 * x_step, -1.0 * y_step),
                            (-1.0 * x_step, y_step),
                            (x_step, y_step)]

        return np.asarray([np.asarray(clockwise_square)])

    # parses a list of M, N-vertex polygon strings, in degrees, centered at the origin, from Treasure Map of the form:
    #   "POLYGON ((ra_1 dec_1, ra_2 dec_2, ... ra_N dec_N))"
    # returns a collection of M, N-vertex polygons, in degrees, centered at the origin of the form:
    #    [
    #      [(ra_1, dec_1), (ra_2, dec_2), ... (ra_N, dec_N)], ...
    #    ]
    @staticmethod
    def get_detector_vertices_from_treasuremap_footprint_string(string_collection):

        output_polygons = []

        for i, s in enumerate(string_collection):
            polygon_vertices = []

            vertex_str = s.split("((")[1].split("))")[0].split(",")
            vertex_tokens = [v.strip().split(" ") for v in vertex_str]

            for v in vertex_tokens:
                polygon_vertices.append((float(v[0]), float(v[1])))

            output_polygons.append(np.asarray(polygon_vertices))

        return np.asarray(output_polygons)

    @staticmethod
    def get_detector_vertices_from_teglon_db(mp_string):

        output_poly = []
        string_list_of_poly = mp_string.replace('MULTIPOLYGON(', '')[:-1].split(")),((")
        for subpoly in string_list_of_poly:
            polygon_vertices = []
            vertex_str = subpoly.replace("((", "").replace("))", "").split(",")
            vertex_tokens = [v.strip().split(" ") for v in vertex_str]

            for v in vertex_tokens:
                polygon_vertices.append((float(v[0]), float(v[1])))
            output_poly.append(np.asarray(polygon_vertices))

        return np.asarray(output_poly)

    def __init__(self,
                 detector_name,
                 detector_vertex_list_collection,
                 detector_width_deg=None,
                 detector_height_deg=None,
                 detector_radius_deg=None,
                 detector_id=None):

        self.id = detector_id
        self.name = detector_name

        # Width, Height, and Radius are holdovers before we started using polygons to represent Detectors directly.
        # They can be included so that they are visible in the database, but radius and area calculation will happen
        # directly from the polygon itself
        self.deg_width = detector_width_deg
        self.deg_height = detector_height_deg
        self.deg_radius = detector_radius_deg
        self.detector_vertex_list_collection = detector_vertex_list_collection

        self.__multipolygon = None
        self.__query_polygon = None
        self.__query_polygon_string = None
        self.__radius_proxy = None

        running_area = 0.0 # sq deg
        for p in self.multipolygon:
            running_area += p.area
        self.area = running_area
        self.__radius_proxy = np.sqrt(self.area/np.pi)


    @property
    def radius_proxy(self):
        return self.__radius_proxy

    @property
    def multipolygon(self):
        if not self.__multipolygon:
            self.__multipolygon = []
            for vertex_list in self.detector_vertex_list_collection:
                self.__multipolygon.append(geometry.Polygon(vertex_list))
        return self.__multipolygon

    @property
    def projected_multipolygon(self):
        return self.multipolygon

    @property
    def query_polygon(self):
        return self.multipolygon

    @property
    def query_polygon_string(self):
        if not self.__query_polygon_string:
            mp_str = "MULTIPOLYGON("
            multipolygon = []
            for geom in self.multipolygon:

                mp = "(("
                # ra_deg, dec_deg = zip(*[(coord_deg[0], coord_deg[1]) for coord_deg in geom.exterior.coords])
                ra_deg, dec_deg = Teglon_Shape.get_coord_lists(geom, convert_radian=False)

                # Here we don't care about emulating lat/lon. These polygons will only be used downstream
                # to instantiate a Tile and project it to a celestial position
                for i in range(len(ra_deg)):
                    mp += "%s %s," % (ra_deg[i], dec_deg[i])

                mp = mp[:-1]  # trim the last ","
                mp += ")),"
                multipolygon.append(mp)

            # Use the multipolygon string to create the WHERE clause
            multipolygon[-1] = multipolygon[-1][:-1]  # trim the last "," from the last object

            for mp in multipolygon:
                mp_str += mp
            mp_str += ")"

            self.__query_polygon_string = mp_str

        return self.__query_polygon_string

    def __str__(self):
        return str(self.__dict__)

