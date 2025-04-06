import 'package:flutter/material.dart';

// Singleton state manager for dashboard
class DashboardState {
  static final DashboardState _instance = DashboardState._internal();
  factory DashboardState() => _instance;
  DashboardState._internal();

  final List<Map<String, dynamic>> completedMolecules = [];
  final List<Map<String, dynamic>> publishedMolecules = [];
  final Set<String> publishedMoleculeNames = {};
  int moleculeCounter = 1;

  void addMolecule(String smiles) {
    completedMolecules.add({
      'name': 'Molecule A$moleculeCounter',
      'smiles': smiles,
      'date': DateTime.now().toString().split(' ')[0],
      'isPublished': false,
    });
    moleculeCounter++;
  }

  void publishMolecule(int index) {
    if (index >= 0 && index < completedMolecules.length) {
      final molecule = completedMolecules[index];

      // Add to published list with published status
      publishedMolecules.add({
        ...molecule,
        'publishedDate': DateTime.now().toString().split(' ')[0],
        'isPublished': true,
      });

      // Track that this molecule has been published
      publishedMoleculeNames.add(molecule['name']);
    }
  }

  bool isMoleculePublished(String name) {
    return publishedMoleculeNames.contains(name);
  }
}

class DashboardScreen extends StatefulWidget {
  const DashboardScreen({Key? key}) : super(key: key);

  @override
  State<DashboardScreen> createState() => _DashboardScreenState();
}

class _DashboardScreenState extends State<DashboardScreen> {
  final dashboardState = DashboardState();
  int selectedTab = 0;

  void addCompletedMolecule(String smiles) {
    setState(() {
      dashboardState.addMolecule(smiles);
    });
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      backgroundColor: Colors.grey[100],
      body: Row(
        children: [
          // Sidebar with shadow
          Card(
            margin: EdgeInsets.zero,
            elevation: 6,
            shape: RoundedRectangleBorder(
              borderRadius: BorderRadius.only(
                topRight: Radius.circular(0),
                bottomRight: Radius.circular(0),
              ),
            ),
            child: Container(
              width: 250,
              color: Colors.white,
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  // Sidebar Header
                  Container(
                    padding: EdgeInsets.symmetric(horizontal: 16, vertical: 16),
                    child: Row(
                      children: [
                        SizedBox(
                          width: 30,
                          height: 30,
                          child: Image.asset(
                            'assets/logo.png',
                            fit: BoxFit.cover,
                          ),
                        ), // Placeholder for logo
                        SizedBox(width: 10),
                        Text(
                          "Medixir AI",
                          style: TextStyle(
                            fontSize: 18,
                            fontWeight: FontWeight.bold,
                          ),
                        ),
                      ],
                    ),
                  ),
                  Divider(height: 1),

                  // New Project Button (Updated)
                  Padding(
                    padding: EdgeInsets.symmetric(horizontal: 16, vertical: 16),
                    child: ElevatedButton.icon(
                      onPressed: () {
                        Navigator.pushNamed(context, '/datasets');
                      },
                      icon: Icon(Icons.add, size: 24, color: Colors.white),
                      label: Text(
                        "New Project",
                        style: TextStyle(fontWeight: FontWeight.bold),
                      ),
                      style: ElevatedButton.styleFrom(
                        backgroundColor: Colors.blue,
                        foregroundColor: Colors.white,
                        padding: EdgeInsets.symmetric(
                          horizontal: 16,
                          vertical: 12,
                        ),
                        minimumSize: Size(double.infinity, 45),
                        shape: RoundedRectangleBorder(
                          borderRadius: BorderRadius.circular(8),
                        ),
                      ),
                    ),
                  ),

                  // Sidebar Navigation Items
                  _sidebarItem(Icons.account_circle, "Account"),
                  _sidebarItem(Icons.people, "Collaborate"),
                  _sidebarItem(Icons.chat, "Chat"),
                  _sidebarItem(Icons.article, "Research Papers"),

                  Spacer(),

                  // Settings at the bottom
                  Divider(height: 1),
                  _sidebarItem(Icons.settings, "Settings"),
                ],
              ),
            ),
          ),

          // Main Content
          Expanded(
            child: Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children: [
                // Top Bar with improved search
                Container(
                  padding: EdgeInsets.all(16),
                  decoration: BoxDecoration(
                    color: Colors.white,
                    boxShadow: [
                      BoxShadow(
                        color: Colors.black.withOpacity(0.05),
                        blurRadius: 5,
                        offset: Offset(0, 2),
                      ),
                    ],
                  ),
                  child: Row(
                    children: [
                      // Improved Search Bar
                      Expanded(
                        child: Container(
                          height: 45,
                          decoration: BoxDecoration(
                            color: Colors.grey[100],
                            borderRadius: BorderRadius.circular(8),
                            border: Border.all(color: Colors.grey[300]!),
                          ),
                          padding: EdgeInsets.symmetric(horizontal: 12),
                          child: Row(
                            children: [
                              Icon(Icons.search, color: Colors.grey[600]),
                              SizedBox(width: 8),
                              Expanded(
                                child: TextField(
                                  decoration: InputDecoration(
                                    hintText: "Search for molecules...",
                                    border: InputBorder.none,
                                    hintStyle: TextStyle(
                                      color: Colors.grey[500],
                                    ),
                                  ),
                                ),
                              ),
                            ],
                          ),
                        ),
                      ),

                      SizedBox(width: 20),

                      // User info and actions
                      Row(
                        children: [
                          Icon(
                            Icons.notifications_none,
                            color: Colors.grey[700],
                          ),
                          SizedBox(width: 16),
                          Text(
                            "User Name",
                            style: TextStyle(
                              fontSize: 16,
                              fontWeight: FontWeight.w500,
                            ),
                          ),
                          SizedBox(width: 12),
                          Container(
                            width: 40,
                            height: 40,
                            decoration: BoxDecoration(
                              color: Colors.grey[300],
                              shape: BoxShape.circle,
                            ),
                            child: Center(
                              child: Text(
                                "U",
                                style: TextStyle(fontWeight: FontWeight.bold),
                              ),
                            ),
                          ), // Placeholder for avatar
                        ],
                      ),
                    ],
                  ),
                ),

                // Main Content Area
                Expanded(
                  child: SingleChildScrollView(
                    padding: EdgeInsets.all(24),
                    child: Column(
                      crossAxisAlignment: CrossAxisAlignment.start,
                      children: [
                        // "Molecules Under Development" Section
                        Text(
                          "Molecules Under Development",
                          style: TextStyle(
                            fontSize: 22,
                            fontWeight: FontWeight.bold,
                          ),
                        ),
                        SizedBox(height: 16),

                        // Molecule Cards in a proper grid layout
                        GridView.count(
                          crossAxisCount: 3,
                          crossAxisSpacing: 20,
                          mainAxisSpacing: 20,
                          shrinkWrap: true,
                          physics: NeverScrollableScrollPhysics(),
                          children: [_moleculeCard("Molecule A-123")],
                        ),

                        SizedBox(height: 32),

                        // Tabs with proper styling
                        Container(
                          decoration: BoxDecoration(
                            border: Border(
                              bottom: BorderSide(color: Colors.grey[300]!),
                            ),
                          ),
                          child: SingleChildScrollView(
                            scrollDirection: Axis.horizontal,
                            child: Row(
                              children: [
                                _tabItem("Completed Molecules", 0),
                                SizedBox(width: 24),
                                _tabItem("Published Molecules", 1),
                              ],
                            ),
                          ),
                        ),

                        SizedBox(height: 20),

                        // List of Completed/Published Molecules
                        selectedTab == 0
                            ? _completedMoleculesList()
                            : _publishedMoleculesList(),
                      ],
                    ),
                  ),
                ),
              ],
            ),
          ),
        ],
      ),
    );
  }

  // Sidebar Item
  Widget _sidebarItem(IconData icon, String title) {
    return InkWell(
      onTap: () {},
      child: Container(
        width: double.infinity,
        padding: EdgeInsets.symmetric(vertical: 16, horizontal: 16),
        child: Row(
          mainAxisSize: MainAxisSize.min,
          children: [
            Icon(icon, color: Colors.grey[800], size: 22),
            SizedBox(width: 16),
            Text(
              title,
              style: TextStyle(fontSize: 16, color: Colors.grey[800]),
            ),
          ],
        ),
      ),
    );
  }

  // Molecule Card with placeholder
  Widget _moleculeCard(String name) {
    return Container(
      decoration: BoxDecoration(
        borderRadius: BorderRadius.circular(8),
        color: Colors.white,
        boxShadow: [
          BoxShadow(
            color: Colors.black.withOpacity(0.05),
            blurRadius: 8,
            offset: Offset(0, 2),
          ),
        ],
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Placeholder for molecule image
          Container(
            height: 150,
            decoration: BoxDecoration(
              color: Colors.grey[200],
              borderRadius: BorderRadius.only(
                topLeft: Radius.circular(8),
                topRight: Radius.circular(8),
              ),
            ),
            child: Center(
              child: Text(
                "Molecule Image",
                style: TextStyle(color: Colors.grey[600]),
              ),
            ),
          ),

          Padding(
            padding: EdgeInsets.all(16),
            child: Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children: [
                Text(
                  name,
                  style: TextStyle(fontSize: 16, fontWeight: FontWeight.bold),
                ),
                SizedBox(height: 8),
                Container(
                  padding: EdgeInsets.symmetric(horizontal: 8, vertical: 4),
                  decoration: BoxDecoration(
                    color: Colors.orange[100],
                    borderRadius: BorderRadius.circular(4),
                    border: Border.all(color: Colors.orange[300]!),
                  ),
                  child: Text(
                    "In Progress",
                    style: TextStyle(
                      fontSize: 12,
                      color: Colors.orange[800],
                      fontWeight: FontWeight.w500,
                    ),
                  ),
                ),
                SizedBox(height: 12),
                Row(
                  mainAxisAlignment: MainAxisAlignment.spaceBetween,
                  children: [
                    Text(
                      "View Details",
                      style: TextStyle(
                        color: Colors.blue,
                        fontWeight: FontWeight.w500,
                      ),
                    ),
                    Icon(Icons.arrow_forward, size: 16, color: Colors.blue),
                  ],
                ),
              ],
            ),
          ),
        ],
      ),
    );
  }

  // Tabs
  Widget _tabItem(String text, int index) {
    final bool isSelected = selectedTab == index;
    return GestureDetector(
      onTap: () => setState(() => selectedTab = index),
      child: Container(
        padding: EdgeInsets.only(bottom: 12),
        decoration: BoxDecoration(
          border: Border(
            bottom: BorderSide(
              color: isSelected ? Colors.blue : Colors.transparent,
              width: 2,
            ),
          ),
        ),
        child: Text(
          text,
          style: TextStyle(
            fontSize: 16,
            fontWeight: isSelected ? FontWeight.bold : FontWeight.w500,
            color: isSelected ? Colors.blue : Colors.grey[700],
          ),
        ),
      ),
    );
  }

  // Completed Molecules List
  Widget _completedMoleculesList() {
    if (dashboardState.completedMolecules.isEmpty) {
      return Center(
        child: Padding(
          padding: EdgeInsets.symmetric(vertical: 40),
          child: Column(
            children: [
              Icon(Icons.science_outlined, size: 48, color: Colors.grey[400]),
              SizedBox(height: 16),
              Text(
                "No molecules completed yet",
                style: TextStyle(
                  fontSize: 16,
                  color: Colors.grey[600],
                  fontWeight: FontWeight.w500,
                ),
              ),
            ],
          ),
        ),
      );
    }

    return Column(
      children:
          dashboardState.completedMolecules.asMap().entries.map((entry) {
            final index = entry.key;
            final molecule = entry.value;
            return _moleculeItem(
              molecule["name"]!,
              molecule["date"]!,
              molecule["smiles"]!,
              index,
            );
          }).toList(),
    );
  }

  // Completed Molecule Item
  Widget _moleculeItem(
    String name,
    String date,
    String smiles, [
    int index = -1,
    bool isPublished = false,
  ]) {
    final bool hasBeenPublished = dashboardState.isMoleculePublished(name);

    return Container(
      margin: EdgeInsets.only(bottom: 12),
      decoration: BoxDecoration(
        borderRadius: BorderRadius.circular(8),
        color: Colors.white,
        border: Border.all(color: Colors.grey[200]!),
      ),
      child: ListTile(
        contentPadding: EdgeInsets.symmetric(horizontal: 16, vertical: 4),
        leading: Container(
          padding: EdgeInsets.all(8),
          decoration: BoxDecoration(
            color: Colors.blue.withOpacity(0.1),
            shape: BoxShape.circle,
          ),
          child: Icon(Icons.science, color: Colors.blue),
        ),
        title: Row(
          children: [
            Text(name, style: TextStyle(fontWeight: FontWeight.bold)),
            SizedBox(width: 8),
            Expanded(
              child: Text(
                smiles,
                style: TextStyle(fontSize: 12, color: Colors.grey[600]),
                overflow: TextOverflow.ellipsis,
              ),
            ),
          ],
        ),
        subtitle: Text("Generated on: $date"),
        trailing: Row(
          mainAxisSize: MainAxisSize.min,
          children: [
            Container(
              padding: EdgeInsets.symmetric(horizontal: 12, vertical: 6),
              decoration: BoxDecoration(
                color: isPublished ? Colors.orange[100] : Colors.green[100],
                borderRadius: BorderRadius.circular(16),
                border: Border.all(
                  color: isPublished ? Colors.orange[300]! : Colors.green[300]!,
                ),
              ),
              child: Text(
                isPublished ? "Published" : "Completed",
                style: TextStyle(
                  color: isPublished ? Colors.orange[800] : Colors.green[800],
                  fontWeight: FontWeight.w500,
                  fontSize: 12,
                ),
              ),
            ),
            // Only show publish button if molecule hasn't been published yet
            if (!isPublished && !hasBeenPublished) ...[
              SizedBox(width: 8),
              ElevatedButton.icon(
                onPressed: () {
                  setState(() {
                    dashboardState.publishMolecule(index);
                  });
                },
                icon: Icon(
                  Icons.publish_rounded,
                  size: 16,
                  color: Colors.white,
                ),
                label: Text(
                  "Publish",
                  style: TextStyle(
                    fontSize: 13,
                    fontWeight: FontWeight.w600,
                    color: Colors.white,
                    letterSpacing: 0.3,
                  ),
                ),
                style: ElevatedButton.styleFrom(
                  backgroundColor: Colors.deepPurple[600],
                  foregroundColor: Colors.white,
                  elevation: 2,
                  shadowColor: Colors.deepPurple.withOpacity(0.3),
                  padding: EdgeInsets.symmetric(horizontal: 16, vertical: 10),
                  minimumSize: Size(0, 0),
                  shape: RoundedRectangleBorder(
                    borderRadius: BorderRadius.circular(20),
                  ),
                ),
              ),
            ],
          ],
        ),
      ),
    );
  }

  // Published Molecules List
  Widget _publishedMoleculesList() {
    if (dashboardState.publishedMolecules.isEmpty) {
      return Center(
        child: Padding(
          padding: EdgeInsets.symmetric(vertical: 40),
          child: Column(
            children: [
              Icon(Icons.science_outlined, size: 48, color: Colors.grey[400]),
              SizedBox(height: 16),
              Text(
                "No molecules published yet",
                style: TextStyle(
                  fontSize: 16,
                  color: Colors.grey[600],
                  fontWeight: FontWeight.w500,
                ),
              ),
            ],
          ),
        ),
      );
    }

    return Column(
      children:
          dashboardState.publishedMolecules.map((molecule) {
            return _moleculeItem(
              molecule["name"]!,
              molecule["date"]!,
              molecule["smiles"]!,
              -1,
              true, // isPublished flag
            );
          }).toList(),
    );
  }
}
