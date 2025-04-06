import 'package:flutter/material.dart';
import 'package:file_picker/file_picker.dart';
import 'package:flutter_dropzone/flutter_dropzone.dart';
import '../utils/api_service.dart';
import 'dart:typed_data';

class DatasetsScreen extends StatefulWidget {
  const DatasetsScreen({super.key});

  @override
  _DatasetsScreenState createState() => _DatasetsScreenState();
}

class _DatasetsScreenState extends State<DatasetsScreen> {
  String? selectedModel; // Stores the selected model
  String? selectedDataset; // Stores the selected dataset
  List<Map<String, String>> uploadedFiles = []; // List of uploaded files
  bool isUploading = false; // Tracks upload state
  bool isProcessing = false; // Tracks processing state
  String? error; // Error message
  bool isDragging = false;
  late DropzoneViewController dropzoneController;
  final Color primarywhite = Color.fromARGB(255, 255, 255, 255);
  final ScrollController _scrollController = ScrollController();

  /// Upload file to the backend
  Future<void> _uploadFile(Uint8List bytes, String fileName) async {
    try {
      bool success = await ApiService.uploadDataset(bytes, fileName);
      if (success) {
        setState(() {
          uploadedFiles.add({
            'name': fileName,
            'size': '${(bytes.lengthInBytes / 1024).toStringAsFixed(2)} KB',
          });
        });
        ScaffoldMessenger.of(
          context,
        ).showSnackBar(SnackBar(content: Text('File uploaded successfully')));
      } else {
        ScaffoldMessenger.of(
          context,
        ).showSnackBar(SnackBar(content: Text('Failed to upload file')));
      }
    } catch (e) {
      ScaffoldMessenger.of(
        context,
      ).showSnackBar(SnackBar(content: Text('Error uploading file: $e')));
    }
  }

  /// Pick and upload a file
  Future<void> _pickAndUploadFile() async {
    try {
      FilePickerResult? result = await FilePicker.platform.pickFiles(
        type: FileType.custom,
        allowedExtensions: ['csv'],
      );

      if (result != null && result.files.isNotEmpty) {
        setState(() {
          isUploading = true;
          error = null;
        });

        final file = result.files.first;
        if (file.bytes != null) {
          await _uploadFile(file.bytes!, file.name);
        }
      }
    } catch (e) {
      setState(() {
        error = e.toString();
      });
    } finally {
      setState(() {
        isUploading = false;
      });
    }
  }

  /// Process the dataset with the selected model
  Future<void> _processDataset() async {
    if (selectedModel == null) {
      ScaffoldMessenger.of(
        context,
      ).showSnackBar(SnackBar(content: Text('Please select a model first')));
      return;
    }

    setState(() {
      isProcessing = true;
    });

    try {
      print('Processing with model: $selectedModel'); // Debug print

      final result = await ApiService.generateMolecule(selectedModel!);
      print('Process result: $result'); // Debug print

      if (result != null) {
        // Navigate regardless of model type
        Navigator.pushNamed(
          context,
          '/generated_molecule',
          arguments: {
            'modelType': selectedModel,
            'modelOutputs': result,
            'moleculeImageUrl': result['molecule_image'],
          },
        );
      } else {
        ScaffoldMessenger.of(
          context,
        ).showSnackBar(SnackBar(content: Text('Failed to generate molecule')));
      }
    } catch (e) {
      print('Error processing dataset: $e');
      ScaffoldMessenger.of(
        context,
      ).showSnackBar(SnackBar(content: Text('Error: $e')));
    } finally {
      setState(() {
        isProcessing = false;
      });
    }
  }

  /// Delete a file from the uploaded list
  void _deleteFile(int index) {
    setState(() {
      if (uploadedFiles[index]["name"] == selectedDataset) {
        selectedDataset = null; // Reset selection if deleted
      }
      uploadedFiles.removeAt(index);
    });
  }

  /// Handle Drag & Drop file addition
  Future<void> _handleDrop(dynamic event) async {
    try {
      final name = event.name as String;
      final size = event.size as int;

      if (!name.toLowerCase().endsWith('.csv')) {
        if (mounted) {
          ScaffoldMessenger.of(context).showSnackBar(
            SnackBar(content: Text("Only CSV files are allowed!")),
          );
        }
        return;
      }

      final bytes = await dropzoneController.getFileData(event);
      await _uploadFile(bytes, name);
    } catch (e) {
      if (mounted) {
        ScaffoldMessenger.of(
          context,
        ).showSnackBar(SnackBar(content: Text("Error handling file: $e")));
      }
    }
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      backgroundColor: Colors.white,
      appBar: AppBar(
        backgroundColor: Colors.blueAccent,
        leading: IconButton(
          icon: Icon(Icons.arrow_back, color: Colors.white),
          onPressed: () {
            Navigator.pop(context);
          },
        ),
      ),
      body: Padding(
        padding: EdgeInsets.all(24.0),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            /// Title & Model Selection & Generate Button
            Row(
              mainAxisAlignment: MainAxisAlignment.spaceBetween,
              children: [
                Text(
                  "Upload Datasets",
                  style: TextStyle(fontSize: 24, fontWeight: FontWeight.bold),
                ),
                Row(
                  children: [
                    /// Model Selection Dropdown
                    Container(
                      decoration: BoxDecoration(
                        border: Border.all(color: Colors.blue, width: 2),
                        borderRadius: BorderRadius.circular(12),
                        color: Colors.white,
                      ),
                      child: DropdownButtonHideUnderline(
                        child: DropdownButton<String>(
                          value: selectedModel,
                          hint: Text("Select Model"),
                          dropdownColor: Colors.white,
                          items:
                              ["Normal", "Biomimetic"].map((model) {
                                return DropdownMenuItem(
                                  value:
                                      model
                                          .toLowerCase(), // Send "normal" or "biomimetic"
                                  child: Padding(
                                    padding: EdgeInsets.symmetric(
                                      horizontal: 14,
                                      vertical: 15,
                                    ),
                                    child: Text(
                                      model,
                                      style: TextStyle(
                                        fontSize: 16,
                                        color: Colors.blue,
                                        fontWeight: FontWeight.w200,
                                      ),
                                    ),
                                  ),
                                );
                              }).toList(),
                          onChanged: (value) {
                            setState(() {
                              selectedModel = value!;
                            });
                          },
                          icon: Icon(Icons.arrow_drop_down, color: Colors.blue),
                        ),
                      ),
                    ),
                    SizedBox(width: 16),

                    /// Generate Button
                    ElevatedButton(
                      onPressed:
                          selectedDataset == null || isUploading || isProcessing
                              ? null
                              : () async {
                                await _processDataset();
                              },
                      style: ElevatedButton.styleFrom(
                        backgroundColor: Colors.blue,
                        foregroundColor: Colors.white,
                        padding: EdgeInsets.symmetric(
                          horizontal: 30,
                          vertical: 18,
                        ),
                        shape: RoundedRectangleBorder(
                          borderRadius: BorderRadius.circular(12),
                        ),
                        disabledBackgroundColor: Colors.grey[300],
                      ),
                      child:
                          isProcessing
                              ? SizedBox(
                                width: 20,
                                height: 20,
                                child: CircularProgressIndicator(
                                  strokeWidth: 2,
                                  valueColor: AlwaysStoppedAnimation<Color>(
                                    Colors.white,
                                  ),
                                ),
                              )
                              : Text(
                                "Generate",
                                style: TextStyle(fontSize: 18),
                              ),
                    ),
                  ],
                ),
              ],
            ),
            SizedBox(height: 16),

            /// Upload Section
            Expanded(
              child: Row(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  /// Drag & Drop Area
                  Expanded(
                    child: Stack(
                      children: [
                        DropzoneView(
                          operation: DragOperation.copy,
                          cursor: CursorType.grab,
                          onCreated:
                              (controller) => dropzoneController = controller,
                          onDrop: _handleDrop,
                          onHover: () => setState(() => isDragging = true),
                          onLeave: () => setState(() => isDragging = false),
                          onError:
                              (String? err) => print('DropzoneError: $err'),
                        ),
                        GestureDetector(
                          onTap: _pickAndUploadFile,
                          child: Container(
                            height: 380,
                            padding: EdgeInsets.all(16),
                            decoration: BoxDecoration(
                              border: Border.all(
                                color:
                                    isDragging
                                        ? Colors.blue
                                        : Colors.grey[300]!,
                                width: 2,
                              ),
                              borderRadius: BorderRadius.circular(12),
                              color:
                                  isDragging ? Color(0xFFF3E5F5) : Colors.white,
                            ),
                            child: Center(
                              child: Column(
                                mainAxisAlignment: MainAxisAlignment.center,
                                children: [
                                  Icon(
                                    Icons.cloud_upload,
                                    size: 52,
                                    color:
                                        isDragging ? primarywhite : Colors.grey,
                                  ),
                                  SizedBox(height: 4),
                                  Text(
                                    "Drag and drop your files here",
                                    style: TextStyle(
                                      fontSize: 22,
                                      fontWeight: FontWeight.w500,
                                    ),
                                    textAlign: TextAlign.center,
                                  ),
                                  SizedBox(height: 4),
                                  Text(
                                    "or",
                                    style: TextStyle(
                                      color: Colors.grey,
                                      fontSize: 21,
                                    ),
                                  ),
                                  SizedBox(height: 4),
                                  Text(
                                    "Browse files",
                                    style: TextStyle(
                                      color: Colors.blue,
                                      fontWeight: FontWeight.bold,
                                      fontSize: 20,
                                    ),
                                  ),
                                ],
                              ),
                            ),
                          ),
                        ),
                      ],
                    ),
                  ),
                  SizedBox(width: 16),

                  /// Uploaded Files List with Selection
                  Expanded(
                    child: Container(
                      height: 380,
                      padding: EdgeInsets.all(16),
                      decoration: BoxDecoration(
                        color: Colors.white,
                        borderRadius: BorderRadius.circular(12),
                        border: Border.all(color: Colors.grey[300]!),
                      ),
                      child: Column(
                        crossAxisAlignment: CrossAxisAlignment.start,
                        children: [
                          Text(
                            "Uploaded Files",
                            style: TextStyle(
                              fontSize: 18,
                              fontWeight: FontWeight.bold,
                            ),
                          ),
                          SizedBox(height: 12),
                          Expanded(
                            child:
                                uploadedFiles.isEmpty
                                    ? Center(
                                      child: Text("No files uploaded yet."),
                                    )
                                    : Scrollbar(
                                      controller: _scrollController,
                                      thumbVisibility: true,
                                      child: ListView.builder(
                                        controller: _scrollController,
                                        itemCount: uploadedFiles.length,
                                        itemBuilder: (context, index) {
                                          return ListTile(
                                            leading: Radio<String>(
                                              value:
                                                  uploadedFiles[index]["name"]!,
                                              groupValue: selectedDataset,
                                              activeColor: Colors.green,
                                              onChanged: (value) {
                                                setState(() {
                                                  selectedDataset = value!;
                                                });
                                              },
                                            ),
                                            title: Text(
                                              uploadedFiles[index]["name"]!,
                                            ),
                                            subtitle: Text(
                                              uploadedFiles[index]["size"]!,
                                            ),
                                            trailing: IconButton(
                                              icon: Icon(
                                                Icons.delete,
                                                color: Colors.red,
                                              ),
                                              onPressed:
                                                  () => _deleteFile(index),
                                            ),
                                          );
                                        },
                                      ),
                                    ),
                          ),
                        ],
                      ),
                    ),
                  ),
                ],
              ),
            ),
          ],
        ),
      ),
    );
  }
}
