from conquest2a.conquest import conquest_coordinates_processor, conquest_input
from conquest2a.writers import xsf_writer

conquest_map = conquest_input({1: "O", 2: "Bi", 3: "Mn", 4: "Mn", 5: "Mn", 6: "Mn"})
path = "./data/test_output_input_coords.in"
coordsproc = conquest_coordinates_processor(path, conquest_map)

xsf_writer("./data/test_output_input_coords.xsf", coordsproc)
