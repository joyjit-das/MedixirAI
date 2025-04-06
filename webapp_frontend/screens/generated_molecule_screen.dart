import 'package:flutter/material.dart';
import '../utils/api_service.dart';
import 'molecule_evaluation_screen.dart';
import 'dashboard_screen.dart'; // Import DashboardScreen

class GeneratedMoleculeScreen extends StatefulWidget {
  final String modelType;
  final Map<String, dynamic>? modelOutputs;
  final String? moleculeImageUrl;

  const GeneratedMoleculeScreen({
    Key? key,
    required this.modelType,
    this.modelOutputs,
    this.moleculeImageUrl,
  }) : super(key: key);

  @override
  State<GeneratedMoleculeScreen> createState() =>
      _GeneratedMoleculeScreenState();
}

class _GeneratedMoleculeScreenState extends State<GeneratedMoleculeScreen> {
  final dashboardState = DashboardState();
  bool isSaved = false;

  Future<void> _saveToDashboard() async {
    if (widget.modelOutputs == null ||
        widget.modelOutputs!['generated_smiles'] == null) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('No molecule data available to save')),
      );
      return;
    }

    dashboardState.addMolecule(widget.modelOutputs!['generated_smiles']);
    setState(() {
      isSaved = true;
    });
    ScaffoldMessenger.of(
      context,
    ).showSnackBar(SnackBar(content: Text('Molecule saved to dashboard')));
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: Text('Generated Molecule (${widget.modelType.toUpperCase()})'),
        backgroundColor: Colors.blue,
        elevation: 0,
        actions: [
          ElevatedButton.icon(
            icon: Icon(Icons.science, color: Colors.white),
            label: Text('Evaluate', style: TextStyle(color: Colors.white)),
            style: ElevatedButton.styleFrom(
              backgroundColor: Colors.green,
              elevation: 0,
              padding: EdgeInsets.symmetric(horizontal: 16),
            ),
            onPressed: () {
              if (widget.modelOutputs?['generated_smiles'] != null) {
                Navigator.push(
                  context,
                  MaterialPageRoute(
                    builder:
                        (context) => MoleculeEvaluationScreen(
                          smiles: widget.modelOutputs!['generated_smiles'],
                        ),
                  ),
                );
              } else {
                ScaffoldMessenger.of(context).showSnackBar(
                  SnackBar(
                    content: Text('No molecule available for evaluation'),
                  ),
                );
              }
            },
          ),
          SizedBox(width: 8),
          ElevatedButton.icon(
            icon: Icon(isSaved ? Icons.check : Icons.save, color: Colors.white),
            label: Text(
              isSaved ? 'Saved' : 'Save',
              style: TextStyle(color: Colors.white),
            ),
            style: ElevatedButton.styleFrom(
              backgroundColor: isSaved ? Colors.green : Colors.teal,
              elevation: 0,
              padding: EdgeInsets.symmetric(horizontal: 16),
            ),
            onPressed: isSaved ? null : _saveToDashboard,
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
      body: SingleChildScrollView(
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.stretch,
          children: [
            // Molecule Image Section
            Container(
              decoration: BoxDecoration(
                //gradient: LinearGradient(
                //  begin: Alignment.topCenter,
                // end: Alignment.bottomCenter,
                //colors: [Colors.blue, Colors.blue.shade50],
                // ),
              ),
              child: Column(
                children: [
                  Padding(
                    padding: EdgeInsets.all(16),
                    child: Text(
                      'Generated Molecule',
                      style: TextStyle(
                        fontSize: 24,
                        fontWeight: FontWeight.bold,
                        color: Colors.white,
                      ),
                    ),
                  ),
                  if (widget.moleculeImageUrl != null)
                    Card(
                      margin: EdgeInsets.all(16),
                      elevation: 8,
                      shape: RoundedRectangleBorder(
                        borderRadius: BorderRadius.circular(15),
                      ),
                      child: Padding(
                        padding: EdgeInsets.all(16),
                        child: Image.network(
                          widget.moleculeImageUrl!,
                          height: 250,
                          fit: BoxFit.contain,
                          errorBuilder: (context, error, stackTrace) {
                            print('Error loading molecule image: $error');
                            return Icon(
                              Icons.broken_image,
                              size: 100,
                              color: Colors.grey,
                            );
                          },
                        ),
                      ),
                    ),
                ],
              ),
            ),

            // Model Outputs Section
            if (widget.modelOutputs != null)
              Padding(
                padding: EdgeInsets.all(16),
                child: Column(
                  crossAxisAlignment: CrossAxisAlignment.start,
                  children: [
                    _buildOutputCard(
                      title: 'SMILES Notation',
                      value: widget.modelOutputs!['generated_smiles'] ?? 'N/A',
                      icon: Icons.science,
                      iconColor: Colors.blue,
                    ),
                    SizedBox(height: 16),
                    Row(
                      children: [
                        Expanded(
                          child: _buildScoreCard(
                            title: 'Similarity',
                            score:
                                widget.modelOutputs!['similarity_score'] ?? 0.0,
                            icon: Icons.compare,
                            color: Colors.green,
                          ),
                        ),
                        SizedBox(width: 16),
                        if (widget.modelType == 'normal')
                          Expanded(
                            child: _buildScoreCard(
                              title: 'Discriminator',
                              score:
                                  widget.modelOutputs!['discriminator_score'] ??
                                  0.0,
                              icon: Icons.analytics,
                              color: Colors.orange,
                            ),
                          ),
                      ],
                    ),
                  ],
                ),
              ),
          ],
        ),
      ),
    );
  }

  Widget _buildOutputCard({
    required String title,
    required String value,
    required IconData icon,
    required Color iconColor,
  }) {
    return Card(
      elevation: 4,
      shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(12)),
      child: Padding(
        padding: EdgeInsets.all(16),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Row(
              children: [
                Icon(icon, color: iconColor),
                SizedBox(width: 8),
                Text(
                  title,
                  style: TextStyle(
                    fontSize: 18,
                    fontWeight: FontWeight.bold,
                    color: Colors.grey[800],
                  ),
                ),
              ],
            ),
            SizedBox(height: 8),
            Container(
              width: double.infinity,
              padding: EdgeInsets.all(12),
              decoration: BoxDecoration(
                color: Colors.grey[100],
                borderRadius: BorderRadius.circular(8),
              ),
              child: Text(
                value,
                style: TextStyle(
                  fontSize: 16,
                  fontFamily: 'Courier',
                  color: Colors.grey[800],
                ),
              ),
            ),
          ],
        ),
      ),
    );
  }

  Widget _buildScoreCard({
    required String title,
    required double score,
    required IconData icon,
    required Color color,
  }) {
    return Card(
      elevation: 4,
      shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(12)),
      child: Padding(
        padding: EdgeInsets.all(16),
        child: Column(
          children: [
            Icon(icon, color: color, size: 32),
            SizedBox(height: 8),
            Text(
              title,
              style: TextStyle(
                fontSize: 16,
                fontWeight: FontWeight.bold,
                color: Colors.grey[800],
              ),
            ),
            SizedBox(height: 4),
            Text(
              score.toStringAsFixed(4),
              style: TextStyle(
                fontSize: 20,
                fontWeight: FontWeight.bold,
                color: color,
              ),
            ),
          ],
        ),
      ),
    );
  }
}
