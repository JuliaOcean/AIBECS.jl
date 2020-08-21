module Rivers

using Unitful
using Unitful: m, s, km, yr
import Unitful: uconvert
using Distances: Haversine
using NearestNeighbors: BruteTree, knn
using OceanGrids
import OceanGrids: regrid

# These river discharge volumetric flow rates and locations were taken from the Dai and Trenberth dataset
# located on UCAR's Research Data Archive (RDA)
# (URL: https://rda.ucar.edu/datasets/ds551.0/index.html#!description)
# and were converted to Julia/AIBECS with permission from Aiguo Dai.
#
# This module only contains the annual mean river discharge of the 200 largest rivers, taken specifically from files `runoff-table-A.txt` and `runoff-table2-top50r.txt` of the original dataset.
#
# Then I copied these files here and wrote them up into a Julia format by hand
# because it was the most straightforward approach for me here...

struct River{T}
    name::String
    lon::Float64
    lat::Float64
    VFR::T # volumetric flow rate
end

# Show method for rivers
function Base.show(io::IO, r::River)
    lat, lon, name, v = r.lat, r.lon, r.name, r.VFR
    print(io, name, " ", (lat,lon), " ", round(Int, ustrip(v)) * unit(v))
end
Base.show(io::IO, m::MIME"text/plain", r::River) = show(io, r)


const RIVERS = [
    River("Amazon",                          -55.5,  -2.0, 6642.0km^3/yr),
    River("Congo",                            15.3,  -4.3, 1308.0km^3/yr),
    River("Orinoco",                         -63.6,   8.1, 1129.0km^3/yr),
    River("Changjiang",                      117.6,  30.8,  944.0km^3/yr),
    River("Brahmaputra",                      89.7,  25.2,  628.0km^3/yr),
    River("Mississippi",                     -90.9,  32.3,  610.0km^3/yr),
    River("Yenisey",                          86.5,  67.4,  599.0km^3/yr),
    River("Paraná",                          -60.7, -32.7,  568.0km^3/yr),
    River("Lena",                            127.4,  70.7,  531.0km^3/yr),
    River("Mekong",                          105.8,  15.1,  525.0km^3/yr),
    River("Tocantins",                       -49.7,  -3.8,  511.0km^3/yr),
    River("Tapajos",                         -56.8,  -5.2,  415.0km^3/yr),
    River("Ob",                               66.6,  66.6,  412.0km^3/yr),
    River("Ganges",                           88.1,  24.5,  404.0km^3/yr),
    River("Irrawaddy",                        96.0,  21.9,  393.0km^3/yr),
    River("St Lawrence",                     -74.7,  45.0,  363.0km^3/yr),
    River("Amur",                            137.0,  50.5,  354.0km^3/yr),
    River("Xingu",                           -52.2,  -3.2,  302.0km^3/yr),
    River("Mackenzie",                      -133.7,  67.5,  290.0km^3/yr),
    River("Xijiang",                         111.3,  23.5,  270.0km^3/yr),
    River("Columbia",                       -121.2,  45.6,  252.0km^3/yr),
    River("Magdalena",                       -74.9,  10.2,  231.0km^3/yr),
    River("Uruguay",                         -58.0, -31.4,  228.0km^3/yr),
    River("Yukon",                          -162.9,  61.9,  212.0km^3/yr),
    River("Atrato",                          -76.7,   6.2,  204.0km^3/yr),
    River("Danube",                           28.7,  45.2,  202.0km^3/yr),
    River("Niger",                             3.5,  11.9,  193.0km^3/yr),
    River("Ogooué",                           10.2,  -0.7,  186.0km^3/yr),
    River("Essequibo",                       -58.6,   5.8,  154.0km^3/yr),
    River("Fraser",                         -121.4,  49.4,  144.0km^3/yr),
    River("Pechora",                          52.2,  67.6,  140.0km^3/yr),
    River("Nelson",                          -97.9,  54.8,  126.0km^3/yr),
    River("Khatanga",                        102.5,  72.0,  124.0km^3/yr),
    River("Sepik",                           142.2,  -4.2,  123.0km^3/yr),
    River("Kolyma",                          158.7,  68.7,  118.0km^3/yr),
    River("Zambeze",                          33.6, -16.1,  117.0km^3/yr),
    River("Severnaya Dvina",                  41.9,  64.1,  112.0km^3/yr),
    River("Indus",                            68.3,  25.4,  104.0km^3/yr),
    River("Sanaga",                           10.1,   3.8,   99.0km^3/yr),
    River("Godavari",                         81.8,  16.9,   97.0km^3/yr),
    River("Rajang",                          112.9,   2.0,   93.0km^3/yr),
    River("Sao Francisco",                   -37.0, -10.0,   90.0km^3/yr),
    River("Usumacinta",                      -91.5,  17.4,   89.0km^3/yr),
    River("Maroni",                          -54.5,   5.0,   86.0km^3/yr),
    River("Rhine",                             6.1,  51.8,   75.0km^3/yr),
    River("Purari",                          145.1,  -7.0,   74.0km^3/yr),
    River("Caniapiscau",                     -69.2,  57.4,   73.0km^3/yr),
    River("Mahanadi",                         84.0,  20.8,   73.0km^3/yr),
    River("Sacramento",                     -121.5,  38.6,   69.0km^3/yr),
    River("Jacui",                           -51.5, -30.0,   69.0km^3/yr),
    River("San Juan CO",                     -77.2,   4.2,   65.0km^3/yr),
    River("Bénoué",                           13.4,   9.3,   64.0km^3/yr),
    River("Kuskokwim",                      -158.1,  61.9,   57.0km^3/yr),
    River("Albany",                          -83.9,  51.3,   57.0km^3/yr),
    River("Huai",                            117.4,  32.9,   57.0km^3/yr),
    River("La Grande",                       -78.6,  53.7,   56.0km^3/yr),
    River("Tigris (Dicle)",                   44.4,  33.3,   56.0km^3/yr),
    River("Krishna",                          80.6,  16.5,   55.0km^3/yr),
    River("Ottawa",                          -76.2,  45.5,   55.0km^3/yr),
    River("Chindwin",                         95.7,  26.0,   55.0km^3/yr),
    River("Taz",                              82.3,  66.6,   55.0km^3/yr),
    River("Po",                               11.6,  44.9,   55.0km^3/yr),
    River("Courantyne",                      -57.7,   5.0,   54.0km^3/yr),
    River("Indigirka",                       147.5,  69.6,   54.0km^3/yr),
    River("Rhone",                             4.6,  43.8,   54.0km^3/yr),
    River("Saguenay",                        -71.6,  48.6,   53.0km^3/yr),
    River("Alabama",                         -87.5,  31.5,   51.0km^3/yr),
    River("Stikine",                        -132.1,  56.7,   51.0km^3/yr),
    River("Eastmain",                        -78.1,  52.2,   49.0km^3/yr),
    River("Brahmani",                         85.0,  21.5,   48.0km^3/yr),
    River("Dnepr",                            35.2,  47.9,   47.0km^3/yr),
    River("Huanghe",                         114.9,  35.2,   47.0km^3/yr),
    River("Jamanxim",                        -55.8,  -5.5,   47.0km^3/yr),
    River("Susquehanna",                     -76.9,  40.2,   46.0km^3/yr),
    River("Beijiang",                        112.9,  23.6,   46.0km^3/yr),
    River("Don",                              40.7,  47.5,   45.0km^3/yr),
    River("Susitna",                        -150.5,  61.5,   45.0km^3/yr),
    River("Narmada",                          73.7,  21.9,   44.0km^3/yr),
    River("Santa Cruz",                      -71.9, -50.3,   41.0km^3/yr),
    River("Oyapock (Oia)",                   -51.9,   3.8,   40.0km^3/yr),
    River("George",                          -65.8,  58.2,   40.0km^3/yr),
    River("Rufiji",                           37.9,  -7.8,   40.0km^3/yr),
    River("Nile",                             31.3,  29.7,   40.0km^3/yr),
    River("Nottaway",                        -77.4,  50.1,   39.0km^3/yr),
    River("Volta",                             0.1,   6.2,   37.0km^3/yr),
    River("Willamette",                     -123.0,  44.9,   37.0km^3/yr),
    River("Skeena",                         -128.4,  54.6,   36.0km^3/yr),
    River("Garonne",                           0.2,  44.4,   36.0km^3/yr),
    River("Red",                             -92.4,  31.3,   36.0km^3/yr),
    River("Pahang",                          102.4,   3.4,   35.0km^3/yr),
    River("Moose",                           -81.3,  50.8,   34.0km^3/yr),
    River("Chao Phraya",                     100.1,  15.3,   34.0km^3/yr),
    River("Cuyuni",                          -58.8,   6.4,   34.0km^3/yr),
    River("Vistula (Wisła)",                  18.8,  54.1,   34.0km^3/yr),
    River("Copper",                         -144.5,  61.5,   34.0km^3/yr),
    River("Saint John",                      -66.8,  46.0,   34.0km^3/yr),
    River("Olenek",                          123.7,  71.9,   34.0km^3/yr),
    River("Rupert",                          -76.9,  51.5,   33.0km^3/yr),
    River("Vuoksi",                           28.8,  61.2,   33.0km^3/yr),
    River("Kamchatka",                       161.6,  56.3,   33.0km^3/yr),
    River("Manicouagan",                     -68.3,  49.2,   33.0km^3/yr),
    River("Esmeraldas",                      -79.4,   0.5,   32.0km^3/yr),
    River("Doce",                            -40.1, -19.4,   32.0km^3/yr),
    River("Yana",                            136.1,  70.8,   32.0km^3/yr),
    River("Jari",                            -52.6,  -0.6,   32.0km^3/yr),
    River("Bío Bío",                         -73.1, -36.8,   32.0km^3/yr),
    River("Severn-CA",                       -88.3,  55.4,   32.0km^3/yr),
    River("Nushagak",                       -157.5,  59.3,   31.0km^3/yr),
    River("Tsiribihina",                      45.0, -19.7,   31.0km^3/yr),
    River("Baker",                           -72.9, -47.3,   30.0km^3/yr),
    River("Araguari",                        -51.4,   0.7,   30.0km^3/yr),
    River("Pur",                              78.2,  67.0,   29.0km^3/yr),
    River("Nass",                           -129.1,  55.2,   29.0km^3/yr),
    River("Douro (Duero)",                    -7.8,  41.2,   29.0km^3/yr),
    River("Negro-AR",                        -63.7, -40.5,   29.0km^3/yr),
    River("Kouilou",                          12.1,  -4.1,   29.0km^3/yr),
    River("Loire",                            -0.8,  47.4,   28.0km^3/yr),
    River("Puelo",                           -72.2, -41.5,   28.0km^3/yr),
    River("Dongjiang",                       114.3,  23.2,   27.0km^3/yr),
    River("Tista",                            89.5,  25.8,   27.0km^3/yr),
    River("Elbe (Labe)",                      11.8,  53.0,   27.0km^3/yr),
    River("Alsek",                          -138.1,  59.4,   27.0km^3/yr),
    River("Tombigbee",                       -87.9,  32.5,   27.0km^3/yr),
    River("Paraiba do Sul",                  -41.3, -21.8,   26.0km^3/yr),
    River("Parnaiba",                        -42.4,  -3.5,   26.0km^3/yr),
    River("Karun",                            48.7,  31.3,   25.0km^3/yr),
    River("Curua-2",                         -54.5,  -5.6,   25.0km^3/yr),
    River("Mezen",                            45.6,  65.0,   25.0km^3/yr),
    River("Kinabatangan",                    117.6,   5.3,   25.0km^3/yr),
    River("Saint-Maurice",                   -72.7,  46.6,   24.0km^3/yr),
    River("Nyanga",                           10.7,  -2.7,   23.0km^3/yr),
    River("Capim",                           -47.8,  -2.5,   23.0km^3/yr),
    River("Grijalva",                        -93.2,  18.0,   22.0km^3/yr),
    River("Senegal",                         -15.5,  16.5,   22.0km^3/yr),
    River("Papaloapan",                      -96.1,  18.2,   22.0km^3/yr),
    River("Euphrates (F",                     44.3,  32.7,   21.0km^3/yr),
    River("Glomma (Glåma)",                   11.1,  59.6,   21.0km^3/yr),
    River("Apalachicola",                    -84.9,  30.7,   20.0km^3/yr),
    River("Nadym",                            72.7,  65.6,   20.0km^3/yr),
    River("Odra (Oder)",                      14.3,  52.8,   20.0km^3/yr),
    River("Daugava (Západnaya Dviná)",        25.9,  56.5,   19.0km^3/yr),
    River("Hayes, MB",                       -92.8,  56.4,   19.0km^3/yr),
    River("Baleine (Great Whale River)",     -77.0,  55.2,   19.0km^3/yr),
    River("Mangoky",                          43.9, -21.8,   19.0km^3/yr),
    River("Feuilles (Rivière aux Feuilles)", -70.4,  58.6,   19.0km^3/yr),
    River("Melezes (Rivière aux Mélèzes)",   -69.6,  57.7,   19.0km^3/yr),
    River("Churchill",                       -94.6,  58.1,   18.0km^3/yr),
    River("Kelantan",                        102.2,   5.8,   18.0km^3/yr),
    River("Clutha",                          169.7, -46.2,   17.0km^3/yr),
    River("Cross",                             9.3,   5.8,   17.0km^3/yr),
    River("Kemi",                             24.7,  65.9,   17.0km^3/yr),
    River("Göta",                             12.4,  58.4,   17.0km^3/yr),
    River("Han",                             127.0,  37.5,   17.0km^3/yr),
    River("Back",                            -96.5,  66.1,   17.0km^3/yr),
    River("Onega",                            38.5,  63.8,   16.0km^3/yr),
    River("Tapi (Tapti)",                     72.9,  21.3,   16.0km^3/yr),
    River("Anabar",                          114.1,  72.0,   16.0km^3/yr),
    River("Jequitinhonh",                    -39.5, -15.9,   14.0km^3/yr),
    River("Kuban",                            38.2,  45.2,   14.0km^3/yr),
    River("Pánuco",                          -98.6,  22.0,   14.0km^3/yr),
    River("Limpopo",                          33.0, -24.5,   14.0km^3/yr),
    River("Winisk",                          -87.2,  54.5,   13.0km^3/yr),
    River("Ebro",                              0.5,  40.8,   13.0km^3/yr),
    River("Narva",                            28.2,  59.4,   12.0km^3/yr),
    River("Colorado-AR",                    -114.6,  32.7,   12.0km^3/yr),
    River("Seine",                             1.2,  49.3,   11.0km^3/yr),
    River("Coopermine",                     -115.4,  67.7,   11.0km^3/yr),
    River("Dnestr",                           29.5,  46.8,   11.0km^3/yr),
    River("Tagus (Tejo)",                     -8.4,  39.5,   11.0km^3/yr),
    River("Itapecuru",                       -44.4,  -3.6,    9.9km^3/yr),
    River("Burdekin",                        147.2, -19.8,   10.0km^3/yr),
    River("Murray",                          142.8, -34.6,    9.4km^3/yr),
    River("Bolshoy Anyuy",                   161.2,  68.2,    8.4km^3/yr),
    River("Grande de Santiago",             -105.1,  21.8,    8.4km^3/yr),
    River("Bandama",                          -4.8,   5.9,    8.3km^3/yr),
    River("Cauvery",                          78.8,  10.8,    7.7km^3/yr),
    River("Juba",                             42.5,   3.8,    7.5km^3/yr),
    River("Fitzroy-QX",                      150.1, -23.1,    7.4km^3/yr),
    River("Brazos",                          -95.8,  29.6,    7.1km^3/yr),
    River("Comoé (Komoé)",                    -3.7,   6.6,    6.8km^3/yr),
    River("Jaguaribe",                       -38.2,  -5.2,    6.6km^3/yr),
    River("Daly",                            130.7, -13.8,    6.3km^3/yr),
    River("Ouémé",                             2.5,   6.9,    5.4km^3/yr),
    River("Guadiana",                         -7.6,  37.8,    5.2km^3/yr),
    River("Anderson",                       -128.4,  68.6,    5.0km^3/yr),
    River("Victoria",                        130.9, -15.5,    4.9km^3/yr),
    River("Roper",                           134.4, -14.7,    4.8km^3/yr),
    River("Liao",                            123.5,  42.2,    4.6km^3/yr),
    River("Orange (Senqu)",                   17.6, -28.8,    4.6km^3/yr),
    River("Paraguacu",                       -39.0, -12.6,    4.3km^3/yr),
    River("Contas [de]",                     -39.3, -14.3,    3.4km^3/yr),
    River("Guadalquivir",                     -6.0,  37.5,    3.3km^3/yr),
    River("Yaqui",                          -109.5,  29.2,    3.3km^3/yr),
    River("Colorado-TX",                     -96.1,  29.3,    2.6km^3/yr),
    River("Rio Grande",                      -97.4,  25.9,    1.5km^3/yr),
    River("de Grey",                         119.2, -20.3,    1.4km^3/yr),
    River("Gascoyne",                        113.8, -24.8,    0.7km^3/yr),
    River("Avon",                            116.1, -31.8,    0.7km^3/yr),
    River("Fortescue",                       116.2, -21.3,    0.3km^3/yr),
    River("Murchison",                       114.5, -27.9,    0.2km^3/yr)
]

citation() = """
- For the dataset itself: Dai, A. 2017. Dai and Trenberth Global River Flow and Continental Discharge Dataset. Research Data Archive at the National Center for Atmospheric Research, Computational and Information Systems Laboratory. https://doi.org/10.5065/D6V69H1T.
- For the formal publication that describe the dataset: Dai, A., and K. E. Trenberth, 2002: Estimates of freshwater discharge from continents: Latitudinal and seasonal variations. J. Hydrometeorol., 3, 660–687.
- For descriptions of updates to the dataset:
    - Dai, A., T. Qian, K. E. Trenberth, and J. D Milliman, 2009: Changes in continental freshwater discharge from 1948-2004. J. Climate, 22, 2773–2791.
    - Dai, A., 2016: Historical and future changes in streamflow and continental runoff: A review, in Terrestrial Water Cycle and Climate Change: Natural and Human-Induced Impacts, Geophysical Monograph 221. Ed. Qiuhong Tang and Taikan Oki, AGU, John Wiley &amp;amp; Sons, 17–37 (DOI: 10.1002/9781118971772).
"""

uconvert(u, r::River) = River(r.name, r.lon, r.lat, uconvert(u, r.VFR))

function load(unit=m^3/s)
    @info """You are about to use the Dai and Trenberth river discharge dataset.
          If you use it for research, please cite:

          $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "Dai_2017", "Dai_Trenberth_2002", "Dai_etal_2009", and "Dai_2016" keys.)
          """
    return uconvert.(unit, RIVERS)
end

function regrid(R::Vector{River{T}}, grd) where T <: Quantity
    lats = [r.lat for r in R]
    lons = [r.lon for r in R]
    depths = zeros(length(R))
    vs = [r.VFR for r in R]
    return regrid(vs, lats, lons, depths, grd)
end
export regrid

end # module

export Rivers
