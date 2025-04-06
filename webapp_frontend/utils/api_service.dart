import 'dart:convert';
import 'package:http/http.dart' as http;
import 'dart:typed_data';
import 'dart:async';
import 'package:flutter/material.dart';

class ApiService {
  static const String baseUrl = 'http://127.0.0.1:5000/api';

  static Future<bool> uploadDataset(
    Uint8List fileBytes,
    String fileName,
  ) async {
    try {
      print('Starting file upload: $fileName');

      var uri = Uri.parse('$baseUrl/upload');
      var request = http.MultipartRequest('POST', uri);

      // Add file to request
      request.files.add(
        http.MultipartFile.fromBytes('file', fileBytes, filename: fileName),
      );

      // Send request
      var streamedResponse = await request.send().timeout(
        Duration(seconds: 30),
        onTimeout: () {
          throw TimeoutException('Upload request timed out');
        },
      );

      // Get response
      var response = await http.Response.fromStream(streamedResponse);
      print('Upload response status: ${response.statusCode}');
      print('Upload response body: ${response.body}');

      if (response.statusCode == 200) {
        print('File uploaded successfully');
        return true;
      } else {
        print('Upload failed with status: ${response.statusCode}');
        return false;
      }
    } catch (e) {
      print('Error uploading file: $e');
      return false;
    }
  }

  static Future<Map<String, dynamic>?> generateMolecule(
    String modelType,
  ) async {
    try {
      print('Generating molecule with model type: $modelType');

      final url = Uri.parse('http://127.0.0.1:5000/api/process');
      final response = await http.post(
        url,
        headers: {'Content-Type': 'application/json'},
        body: jsonEncode({
          'model': modelType,
          'show_plots': false, // Add this parameter
        }),
      );

      print('Generate Molecule - Response status: ${response.statusCode}');
      print('Generate Molecule - Response body: ${response.body}');

      if (response.statusCode == 200) {
        final data = jsonDecode(response.body);
        print('Decoded response data: $data'); // Debug print

        // Validate minimum required fields
        if (!data.containsKey('molecule_image') ||
            !data.containsKey('generated_smiles')) {
          print('Missing required fields in response');
          return null;
        }

        // Construct full URLs for images
        if (data['molecule_image'] != null) {
          data['molecule_image'] =
              'http://127.0.0.1:5000${data['molecule_image']}';
        }

        // Handle model-specific fields
        if (modelType.toLowerCase() == 'biomimetic') {
          // Ensure biomimetic-specific fields have default values if missing
          data['similarity_score'] ??= 0.0;
          data['discriminator_score'] ??= 0.0;
        }

        print('Processed and validated data: $data'); // Debug print
        return data;
      } else {
        print('Failed to process dataset: ${response.statusCode}');
        return null;
      }
    } catch (e) {
      print('Error in generateMolecule: $e');
      return null;
    }
  }

  // Helper method to validate response
  static bool _validateMoleculeResponse(Map<String, dynamic> data) {
    final requiredFields = [
      'closest_smiles',
      'similarity_score',
      'discriminator_score',
      'vector',
      'molecule_image',
    ];

    for (final field in requiredFields) {
      if (!data.containsKey(field) || data[field] == null) {
        print('Missing required field: $field');
        return false;
      }
    }

    // Validate specific field types
    try {
      String _ = data['closest_smiles'] as String;
      double _ = data['similarity_score'] as double;
      double _ = data['discriminator_score'] as double;
      List _ = data['vector'] as List;
      String _ = data['molecule_image'] as String;
      return true;
    } catch (e) {
      print('Invalid field type in response: $e');
      return false;
    }
  }

  static Future<Map<String, dynamic>?> explainGeneration(
    List<double> vector,
  ) async {
    try {
      final response = await http.post(
        Uri.parse('$baseUrl/explain'),
        headers: {"Content-Type": "application/json"},
        body: jsonEncode({"vector": vector}),
      );
      if (response.statusCode == 200) {
        return jsonDecode(response.body);
      } else {
        print('Error: ${response.body}');
        return null;
      }
    } catch (e) {
      print('Error: $e');
      return null;
    }
  }

  static Future<bool> processDataset(
    String datasetName,
    String modelType,
  ) async {
    try {
      print('Processing dataset: $datasetName with model: $modelType');

      final response = await http.post(
        Uri.parse('$baseUrl/process'),
        headers: {'Content-Type': 'application/json'},
        body: json.encode({'dataset': datasetName, 'model': modelType}),
      );

      print('Process response: ${response.body}');

      if (response.statusCode == 200) {
        final responseData = json.decode(response.body);
        return responseData['status'] == 'success';
      } else {
        print('Failed to process dataset: ${response.body}');
        return false;
      }
    } catch (e) {
      print('Error processing dataset: $e');
      return false;
    }
  }

  static Future<Map<String, dynamic>> evaluateMolecule(String smiles) async {
    try {
      final url = Uri.parse('http://127.0.0.1:5000/api/evaluate');
      final response = await http.post(
        url,
        headers: {'Content-Type': 'application/json'},
        body: jsonEncode({'smiles': smiles}),
      );

      if (response.statusCode == 200) {
        return jsonDecode(response.body);
      } else {
        throw Exception('Failed to evaluate molecule: ${response.statusCode}');
      }
    } catch (e) {
      throw Exception('Error evaluating molecule: $e');
    }
  }
}

class MoleculePage extends StatefulWidget {
  @override
  _MoleculePageState createState() => _MoleculePageState();
}

class _MoleculePageState extends State<MoleculePage> {
  String? shapGraphUrl;
  String? moleculeImageUrl;
  String? error;
  bool isLoading = true;

  Future<void> processDataset() async {
    try {
      setState(() {
        isLoading = true;
        error = null;
      });

      final result = await ApiService.generateMolecule('normal');
      print('Process Dataset - Result: $result');

      if (result != null &&
          result.containsKey('shap_graph') &&
          result.containsKey('molecule_image')) {
        setState(() {
          shapGraphUrl = result['shap_graph'];
          moleculeImageUrl = result['molecule_image'];
          isLoading = false;
        });
      } else {
        setState(() {
          error = 'Invalid response from server';
          isLoading = false;
        });
      }
    } catch (e) {
      print('Error processing dataset: $e');
      setState(() {
        error = 'Error: $e';
        isLoading = false;
      });
    }
  }

  @override
  void initState() {
    super.initState();
    processDataset();
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: Text('Generated Molecule'),
        actions: [
          IconButton(icon: Icon(Icons.refresh), onPressed: processDataset),
        ],
      ),
      body:
          isLoading
              ? Center(child: CircularProgressIndicator())
              : error != null
              ? Center(
                child: Column(
                  mainAxisAlignment: MainAxisAlignment.center,
                  children: [
                    Text(
                      error!,
                      style: TextStyle(color: Colors.red),
                      textAlign: TextAlign.center,
                    ),
                    SizedBox(height: 16),
                    ElevatedButton(
                      onPressed: processDataset,
                      child: Text('Retry'),
                    ),
                  ],
                ),
              )
              : SingleChildScrollView(
                padding: EdgeInsets.all(16),
                child: Column(
                  children: [
                    if (shapGraphUrl != null) ...[
                      Text(
                        'SHAP Graph',
                        style: TextStyle(
                          fontSize: 18,
                          fontWeight: FontWeight.bold,
                        ),
                      ),
                      Card(
                        child: Image.network(
                          shapGraphUrl!,
                          errorBuilder: (context, error, stackTrace) {
                            print('Error loading SHAP graph: $error');
                            return Icon(Icons.broken_image, size: 100);
                          },
                        ),
                      ),
                      SizedBox(height: 20),
                    ],
                    if (moleculeImageUrl != null) ...[
                      Text(
                        'Generated Molecule',
                        style: TextStyle(
                          fontSize: 18,
                          fontWeight: FontWeight.bold,
                        ),
                      ),
                      Card(
                        child: Image.network(
                          moleculeImageUrl!,
                          errorBuilder: (context, error, stackTrace) {
                            print('Error loading molecule image: $error');
                            return Icon(Icons.broken_image, size: 100);
                          },
                        ),
                      ),
                    ],
                  ],
                ),
              ),
    );
  }
}
