import Extractor as cyto
import pandas as pd
import click

@click.command()
@click.option('--file', help='The input file. Must be a csv file with columns "Cytogenetics" and "ID".')
@click.option('--version', help='The config version to run')
@click.option('--bool_mode', help='Return type for boolean values. If "string" then convert to strings, otherwise return bools.', default='string')
@click.option('--show_all', help='whether to show all tested abnormalities or only those detected for at least one sample in the output columns', default=True)
@click.option('--fill', help='By default the extractor returns nothing when an abnormality is not found. If fill=True this will be returned as False instead.', default=True)
@click.option('--output', help='output file name. If None it will be [file]_extracted.csv', default=None)
def hello(file, version, bool_mode, show_all, fill, output):
	"""Run the extraction for an input file"""
	if version not in configs:
		click.echo('The specified version is not available')
		return
	
	try:
		df = pd.read_csv(file)
	except:
		click.echo("Input file not found.")
		return
	
	results = [cyto.extract_from_string(x, configs[version], bool_mode=bool_mode) for x in df['Cytogenetics']]
	rows = []
	copy = ['error', 'error_message', 'Warnings', 'fish_available']
	
	for res in results:
		r = {x: res[x] for x in copy}
		for k, v in res['result'].items():
			r[k] = v
		rows.append(r)

	df2 = pd.DataFrame(rows)

	if show_all:
		all_abn = list(configs[version].values())
		missing = [x for x in all_abn if x not in df2.columns]
		for x in missing:
			df2[x] = False
	
	if fill:
		df2 = df2.fillna(False)
	
	df = df.join(df2, rsuffix='_result')

	if output == None:
		output = file[:-4] + '_extracted.csv'

	df.to_csv(output)

configs = cyto.available_configs()

if __name__ == '__main__':
	hello()