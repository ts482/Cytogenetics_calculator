from flask import Flask, request
from flask_cors import CORS
import Extractor as cyto

app = Flask(__name__)
CORS(app)

#setup
abnormalities = cyto.base_extraction()
props = cyto.setup(abnormalities)



@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"

@app.route('/karyotype/extract', methods=['POST'])
@app.route('/extract/karyotype/extract', methods=['post']) #for rosalind proxy
def karyotype_api():
    input = request.json
    if 'fish' not in input:
        input['fish'] = None
    bool_mode = 'string'
    if 'bool_mode' in input:
        bool_mode = input['bool_mode']
    extracted = cyto.extract_from_string(input['karyotype_string'], props, bool_mode=bool_mode, fish = input['fish'])
    return extracted

# TODO bulk API to extract multiple strings at once

# could easily add an endpoint that also allows the abnormalities for extraction to be specified
# rather than using the app-level default 

@app.route('/karyotype/ping', methods=['POST'])
def ping():
    input = request.json
    return input