from flask import Flask, request
from flask_cors import CORS
import Extractor as cyto

DEFAULT_VERSION = "BJH2021"

app = Flask(__name__)
CORS(app)

#setup
configs = cyto.available_configs()

@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"

@app.route('/karyotype/extract', methods=['POST'])
@app.route('/extract/karyotype/extract', methods=['POST']) #for rosalind proxy
def karyotype_api():
    input = request.json
    if 'fish' not in input:
        input['fish'] = None
    bool_mode = 'string'
    if 'bool_mode' in input:
        bool_mode = input['bool_mode']

    version = DEFAULT_VERSION
    if 'version' in input:
        version = input['version']
    if version not in configs:
        return {'error': True, 'error_message': ['model version not recognised']}
    
    extracted = cyto.extract_from_string(input['karyotype_string'], configs[version], bool_mode=bool_mode, fish = input['fish'])
    return extracted

# TODO bulk API to extract multiple strings at once

# could easily add an endpoint that also allows the abnormalities for extraction to be specified
# rather than using the app-level default 

@app.route('/karyotype/ping', methods=['POST'])
def ping():
    input = request.json
    return input