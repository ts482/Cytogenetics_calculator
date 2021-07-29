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
def karyotype_api():
    input = request.json
    extracted = cyto.extract_from_string(input['karyotype_string'], props)
    return extracted

@app.route('/karyotype/ping', methods=['POST'])
def ping():
    input = request.json
    return input