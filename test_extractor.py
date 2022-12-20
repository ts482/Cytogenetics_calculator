import pytest
import Extractor as cyto
import pandas as pd
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

df = pd.read_csv('test_data/karyotypes.csv')
outputs = ['error', 'error_message', 'Warnings', 'fish_available', 'Error', 'Error description', 'Number of cytogenetic abnormalities', 'Monosomy', 'Polysomy', 'Structural', 'inv(3)', '-X', 'Monosomy7', 't(12p)', 't(3;21)', 'del5q', 't(5q)', 't(3;5)', 't(1;3)', 't(3;3)', 'abnormal(17p)', 't(17p)', 'del12p', 't(v;11)', 't(9;11)', 't(6;11)', '-Y', 'del11q', 'del13q', 'del7q', 'idic(X)(q13)', 'isochromosome17q', 'Monosomy13', 'Monosomy17', 'Monosomy5', 't(11;16)(q23.3;p13.3)', 't(2;11)', 't(5;10)', 't(5;12)', 't(5;17)', 't(5;7)', 't(1;22)', 't(6;9)', 't(9;22)', 't(16;16)', 'inv(16)', 't(8;21)', 't(15;17)', 't(10;11)', 't(8;16)(p11;p13)', 't(3q26.2;v)']
results = ['Error', 'Number of cytogenetic abnormalities', 'Monosomy', 'Polysomy', 'Structural', 'inv(3)', '-X', 'Monosomy7', 't(12p)', 't(3;21)', 'del5q', 't(5q)', 't(3;5)', 't(1;3)', 't(3;3)', 'abnormal(17p)', 't(17p)', 'del12p', 't(v;11)', 't(9;11)', 't(6;11)', '-Y', 'del11q', 'del13q', 'del7q', 'idic(X)(q13)', 'isochromosome17q', 'Monosomy13', 'Monosomy17', 'Monosomy5', 't(11;16)(q23.3;p13.3)', 't(2;11)', 't(5;10)', 't(5;12)', 't(5;17)', 't(5;7)', 't(1;22)', 't(6;9)', 't(9;22)', 't(16;16)', 'inv(16)', 't(8;21)', 't(15;17)', 't(10;11)', 't(8;16)(p11;p13)', 't(3q26.2;v)']
@pytest.mark.parametrize("test_input,expected", [(row['Cytogenetics'], row[outputs]) for index, row in df.iterrows()])
def test_file(test_input, expected):
	bool_mode = None
	extracted = cyto.extract_from_string(test_input, configs[test_version], bool_mode=bool_mode, fish = None)
	assert extracted['error'] == expected['error']
	assert extracted['fish_available'] == expected['fish_available']
	for k in results:
		if k in extracted['result']:
			assert extracted['result'][k] == expected[k]
		else:
			assert expected[k] == False
	
	#this is not a great test but ok for now
	assert str(extracted['error_message']) == expected['error_message']
	assert str(extracted['result']['Error description']) == expected['Error description']
	assert str(extracted['result']['Warnings']) == expected['Warnings']