import Extractor as cyto
configs = cyto.available_configs()
test_version = 'ELN2022'

def test_configs():
	configs = cyto.available_configs()
	assert 'BJH2021' in configs
	assert 'ELN2022' in configs

#extracted = cyto.extract_from_string(
# input['karyotype_string'], 
# configs[version], 
# bool_mode=bool_mode, 
# fish = input['fish']
# )

def test_normal():
	karyotype = "46,XY[10]"
	bool_mode = 'string'
	extracted = cyto.extract_from_string(karyotype, configs[test_version], bool_mode=bool_mode, fish = None)
	assert extracted['error'] == False
	assert len(extracted['error_message']) == 0
	assert len(extracted['Warnings']) == 0 
	assert extracted['fish_available'] == False
	assert 'result' in extracted
	assert extracted['result']['Error'] == False
	assert extracted['result']['Number of cytogenetic abnormalities'] == 0
	assert extracted['result']['Monosomy'] == 0
	assert extracted['result']['Polysomy'] == 0
	assert extracted['result']['Structural'] == 0
	assert len(extracted['result']['Error description']) == 0
	assert len(extracted['result']['Warnings']) == 0