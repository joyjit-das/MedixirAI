import 'package:flutter/material.dart';
import 'package:google_fonts/google_fonts.dart';
import '../utils/api_service.dart';

class MoleculeEvaluationScreen extends StatefulWidget {
  final String smiles;

  const MoleculeEvaluationScreen({Key? key, required this.smiles})
    : super(key: key);

  @override
  _MoleculeEvaluationScreenState createState() =>
      _MoleculeEvaluationScreenState();
}

class _MoleculeEvaluationScreenState extends State<MoleculeEvaluationScreen> {
  bool isLoading = true;
  String? error;
  Map<String, dynamic>? metrics;
  String? report;

  @override
  void initState() {
    super.initState();
    _loadEvaluation();
  }

  Future<void> _loadEvaluation() async {
    try {
      setState(() {
        isLoading = true;
        error = null;
      });

      final results = await ApiService.evaluateMolecule(widget.smiles);

      setState(() {
        metrics = results['metrics'];
        report = results['report'];
        isLoading = false;
      });
    } catch (e) {
      setState(() {
        error = e.toString();
        isLoading = false;
      });
    }
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: Text('Molecule Evaluation'),
        backgroundColor: Colors.blue,
        elevation: 0,
        actions: [
          ElevatedButton.icon(
            icon: Icon(Icons.refresh, color: Colors.white),
            label: Text('Retry', style: TextStyle(color: Colors.white)),
            style: ElevatedButton.styleFrom(
              backgroundColor: Colors.green,
              elevation: 0,
              padding: EdgeInsets.symmetric(horizontal: 16),
            ),
            onPressed: _loadEvaluation,
          ),
          SizedBox(width: 8),
          ElevatedButton.icon(
            icon: Icon(Icons.dashboard, color: Colors.white),
            label: Text('Dashboard', style: TextStyle(color: Colors.white)),
            style: ElevatedButton.styleFrom(
              backgroundColor: Colors.orange,
              elevation: 0,
              padding: EdgeInsets.symmetric(horizontal: 16),
            ),
            onPressed:
                () => Navigator.pushReplacementNamed(context, '/dashboard'),
          ),
          SizedBox(width: 16),
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
                    Icon(Icons.error_outline, size: 48, color: Colors.red),
                    SizedBox(height: 16),
                    Text(error!, style: TextStyle(color: Colors.red)),
                    SizedBox(height: 16),
                    ElevatedButton(
                      onPressed: _loadEvaluation,
                      child: Text('Retry'),
                    ),
                  ],
                ),
              )
              : SingleChildScrollView(
                padding: EdgeInsets.all(16),
                child: Column(
                  children: [
                    Card(
                      elevation: 4,
                      shape: RoundedRectangleBorder(
                        borderRadius: BorderRadius.circular(12),
                      ),
                      child: Padding(
                        padding: EdgeInsets.all(24),
                        child: Column(
                          crossAxisAlignment: CrossAxisAlignment.start,
                          children: [
                            Text(
                              'Evaluation Report',
                              style: TextStyle(
                                fontSize: 28,
                                fontWeight: FontWeight.bold,
                                color: Colors.blue[900],
                                letterSpacing: 0.5,
                              ),
                            ),
                            SizedBox(height: 24),
                            _buildReportContent(
                              report ?? "No report available.",
                            ),
                          ],
                        ),
                      ),
                    ),
                  ],
                ),
              ),
    );
  }

  Widget _buildReportContent(String reportText) {
    final sections = reportText.split('\n\n');
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children:
          sections.map((section) {
            if (section.trim().isEmpty) return SizedBox.shrink();

            final lines = section.trim().split('\n');
            if (lines.isEmpty) return SizedBox.shrink();

            // Check if this is a section header
            if (lines[0].toUpperCase() == lines[0] && !lines[0].contains(':')) {
              return Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  SizedBox(height: 24),
                  Text(
                    lines[0],
                    style: GoogleFonts.roboto(
                      fontSize: 20,
                      fontWeight: FontWeight.bold,
                      color: Colors.blue[900],
                      letterSpacing: 0.5,
                    ),
                  ),
                  SizedBox(height: 16),
                  ...lines.skip(1).map((line) => _buildReportLine(line)),
                ],
              );
            }

            return Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children:
                  lines.map((line) {
                    // Handle numeric value checks (LogP, HBonds, etc.)
                    if (line.contains(':')) {
                      final parts = line.split(':');
                      if (parts.length == 2) {
                        // Clean up any trailing [] from numeric values
                        final cleanedValue = parts[1].replaceAll(
                          RegExp(r'\s*\[\s*\]\s*'),
                          '',
                        );
                        return _buildReportLine('${parts[0]}:$cleanedValue');
                      }
                    }

                    // Handle status badges
                    if (line.contains('[') && line.contains(']')) {
                      final match = RegExp(
                        r'(.*?)\[(.*?)\](.*)',
                      ).firstMatch(line);
                      if (match != null) {
                        final beforeStatus = match.group(1)?.trim() ?? '';
                        final status = match.group(2)?.trim() ?? '';
                        if (_isValidStatus(status)) {
                          return _buildReportLine('$beforeStatus[$status]');
                        }
                      }
                    }
                    return _buildReportLine(line);
                  }).toList(),
            );
          }).toList(),
    );
  }

  Widget _buildReportLine(String line) {
    if (line.trim().isEmpty) return SizedBox(height: 8);

    // Check if it's a metric line with a status
    if (line.contains('[') && line.contains(']')) {
      final match = RegExp(r'\[(.*?)\]').firstMatch(line);
      if (match != null) {
        final status = match.group(1)?.trim() ?? '';
        final beforeStatus = line.substring(0, line.indexOf('[')).trim();

        // Only show status if it's a valid one
        if (_isValidStatus(status)) {
          return Padding(
            padding: EdgeInsets.symmetric(vertical: 4),
            child: Row(
              children: [
                Expanded(
                  child: Text(
                    beforeStatus,
                    style: GoogleFonts.roboto(fontSize: 16, height: 1.5),
                  ),
                ),
                Container(
                  padding: EdgeInsets.symmetric(horizontal: 12, vertical: 4),
                  decoration: BoxDecoration(
                    color: _getStatusColor(status),
                    borderRadius: BorderRadius.circular(4),
                  ),
                  child: Text(
                    status, // Removed the brackets here
                    style: GoogleFonts.roboto(
                      color: Colors.white,
                      fontWeight: FontWeight.bold,
                      fontSize: 14,
                    ),
                  ),
                ),
              ],
            ),
          );
        }
      }
    }

    // Regular line
    return Padding(
      padding: EdgeInsets.symmetric(vertical: 4),
      child: Text(line, style: GoogleFonts.roboto(fontSize: 16, height: 1.5)),
    );
  }

  bool _isValidStatus(String status) {
    final validStatuses = {
      'PASS',
      'FAIL',
      'LOW RISK',
      'MEDIUM RISK',
      'HIGH RISK',
      'DRUG-LIKE',
      'NON-DRUG-LIKE',
    };
    return validStatuses.contains(status.toUpperCase().trim());
  }

  Color _getStatusColor(String status) {
    switch (status.toUpperCase().trim()) {
      case 'PASS':
      case 'DRUG-LIKE':
      case 'LOW RISK':
        return Colors.green[700]!;

      case 'FAIL':
      case 'NON-DRUG-LIKE':
      case 'HIGH RISK':
        return Colors.red[700]!;

      case 'MEDIUM RISK':
        return Colors.orange[700]!;

      default:
        return Colors
            .blue[700]!; // This should never happen due to _isValidStatus check
    }
  }
}
