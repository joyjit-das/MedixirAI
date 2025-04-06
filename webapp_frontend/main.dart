import 'package:flutter/material.dart';
import 'package:firebase_core/firebase_core.dart';
import 'firebase_options.dart';
import 'package:firebase_auth/firebase_auth.dart';
import 'screens/login_screen.dart';
import 'screens/dashboard_screen.dart';
import 'screens/datasets_screen.dart';
import 'screens/generated_molecule_screen.dart';
import 'screens/molecule_evaluation_screen.dart';
import 'package:google_fonts/google_fonts.dart';

void main() async {
  WidgetsFlutterBinding.ensureInitialized();
  await Firebase.initializeApp(options: DefaultFirebaseOptions.currentPlatform);
  runApp(MyApp());
}

class MyApp extends StatelessWidget {
  const MyApp({super.key});

  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      debugShowCheckedModeBanner: false,
      title: 'MedixirAI',
      theme: ThemeData(
        primarySwatch: Colors.blue,
        visualDensity: VisualDensity.adaptivePlatformDensity,
        textTheme: GoogleFonts.sourceCodeProTextTheme(
          Theme.of(context).textTheme,
        ),
      ),
      initialRoute: '/',
      routes: {
        '/': (context) => LoginScreen(),
        '/dashboard': (context) => DashboardScreen(),
        '/datasets': (context) => DatasetsScreen(),
        '/generated_molecule': (context) {
          final args =
              ModalRoute.of(context)?.settings.arguments
                  as Map<String, dynamic>?;
          return GeneratedMoleculeScreen(
            modelType: args?['modelType'] ?? 'normal',
            modelOutputs: args?['modelOutputs'],
            moleculeImageUrl: args?['moleculeImageUrl'], // Remove shapGraphUrl
          );
        },
        '/molecule_evaluation':
            (context) => MoleculeEvaluationScreen(smiles: ''),
      },
      // Handle navigation when user is not authenticated
      onGenerateRoute: (settings) {
        // If user tries to access any route other than login while not authenticated,
        // redirect to login screen
        if (!isAuthenticated() && settings.name != '/') {
          return MaterialPageRoute(builder: (context) => LoginScreen());
        }
        return null;
      },
    );
  }
}

// Helper function to check authentication status
bool isAuthenticated() {
  final auth = FirebaseAuth.instance;
  return auth.currentUser != null;
}
