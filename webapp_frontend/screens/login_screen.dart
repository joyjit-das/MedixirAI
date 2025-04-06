import 'package:flutter/material.dart';
import 'dashboard_screen.dart';

class LoginScreen extends StatefulWidget {
  const LoginScreen({super.key});

  @override
  _LoginScreenState createState() => _LoginScreenState();
}

class _LoginScreenState extends State<LoginScreen> {
  final TextEditingController _emailController = TextEditingController();
  final TextEditingController _passwordController = TextEditingController();

  String selectedRole = ""; // Holds selected role ("Researcher" or "Admin")

  bool isValidEmail(String email) {
    return RegExp(
      r"^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$",
    ).hasMatch(email);
  }

  void showCustomSnackbar(BuildContext context, String message, Color bgColor) {
    ScaffoldMessenger.of(context).showSnackBar(
      SnackBar(
        content: Text(
          message,
          style: TextStyle(
            fontSize: 16,
            fontWeight: FontWeight.bold,
            color: Colors.white,
          ),
        ),
        backgroundColor: bgColor,
        behavior: SnackBarBehavior.floating,
        margin: EdgeInsets.all(20),
      ),
    );
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      backgroundColor: Colors.white,
      body: Center(
        child: Container(
          width: 400,
          padding: EdgeInsets.all(20),
          decoration: BoxDecoration(
            color: Colors.white,
            borderRadius: BorderRadius.circular(12),
            boxShadow: [
              BoxShadow(color: Colors.black12, blurRadius: 10, spreadRadius: 2),
            ],
          ),
          child: Column(
            mainAxisSize: MainAxisSize.min,
            children: [
              Image.asset('assets/logo.png', width: 100, height: 100),
              SizedBox(height: 5),
              Text(
                "Medixir AI",
                style: TextStyle(fontSize: 22, fontWeight: FontWeight.bold),
              ),
              SizedBox(height: 20),

              // Role Selection
              Row(
                mainAxisAlignment: MainAxisAlignment.spaceBetween,
                children: [
                  Expanded(
                    child: ElevatedButton.icon(
                      onPressed:
                          () => setState(() => selectedRole = "Researcher"),
                      icon: Icon(Icons.person, color: Colors.white),
                      label: Text(
                        "Researcher",
                        style: TextStyle(
                          fontSize: 17,
                          fontWeight: FontWeight.bold,
                          color: Colors.white,
                        ),
                      ),
                      style: ElevatedButton.styleFrom(
                        backgroundColor:
                            selectedRole == "Researcher"
                                ? Colors.deepPurple
                                : Colors.green,
                        padding: EdgeInsets.symmetric(vertical: 16),
                      ),
                    ),
                  ),
                  SizedBox(width: 10),
                  Expanded(
                    child: ElevatedButton.icon(
                      onPressed: () => setState(() => selectedRole = "Admin"),
                      icon: Icon(
                        Icons.admin_panel_settings,
                        color: Colors.white,
                      ),
                      label: Text(
                        "Admin",
                        style: TextStyle(
                          fontSize: 17,
                          fontWeight: FontWeight.bold,
                          color: Colors.white,
                        ),
                      ),
                      style: ElevatedButton.styleFrom(
                        backgroundColor:
                            selectedRole == "Admin"
                                ? Colors.deepPurple
                                : Colors.green,
                        padding: EdgeInsets.symmetric(vertical: 16),
                      ),
                    ),
                  ),
                ],
              ),
              SizedBox(height: 15),

              TextField(
                controller: _emailController,
                decoration: InputDecoration(
                  labelText: "Email",
                  border: OutlineInputBorder(),
                ),
              ),
              SizedBox(height: 10),
              TextField(
                controller: _passwordController,
                obscureText: true,
                decoration: InputDecoration(
                  labelText: "Password",
                  border: OutlineInputBorder(),
                ),
              ),
              SizedBox(height: 10),

              Align(
                alignment: Alignment.centerRight,
                child: TextButton(
                  onPressed: () {}, // Implement forgot password logic
                  child: Text(
                    "Forgot password?",
                    style: TextStyle(fontSize: 15, color: Colors.deepPurple),
                  ),
                ),
              ),
              SizedBox(height: 10),

              // Sign in button (Manual Login)
              ElevatedButton.icon(
                onPressed: () {
                  if (selectedRole.isEmpty) {
                    showCustomSnackbar(
                      context,
                      "Please select Researcher or Admin.",
                      Colors.redAccent,
                    );
                  } else {
                    Navigator.pushReplacementNamed(context, '/dashboard');
                  }
                },
                icon: Icon(Icons.login, color: Colors.white),
                label: Text(
                  "Sign in",
                  style: TextStyle(
                    fontSize: 17,
                    fontWeight: FontWeight.bold,
                    color: Colors.white,
                  ),
                ),
                style: ElevatedButton.styleFrom(
                  backgroundColor: Colors.deepPurple,
                  padding: EdgeInsets.symmetric(vertical: 16, horizontal: 30),
                ),
              ),
              SizedBox(height: 10),

              // Google Sign-In Button
              Column(
                children: [
                  Text(
                    "or",
                    style: TextStyle(
                      fontSize: 16,
                      fontWeight: FontWeight.bold,
                      color: Colors.grey,
                    ),
                  ),
                  SizedBox(height: 10),
                  ElevatedButton(
                    onPressed: () {
                      Navigator.push(
                        context,
                        MaterialPageRoute(
                          builder: (context) => DashboardScreen(),
                        ),
                      );
                    },
                    style: ElevatedButton.styleFrom(
                      backgroundColor: Colors.white,
                      padding: EdgeInsets.symmetric(
                        vertical: 16,
                        horizontal: 30,
                      ),
                      side: BorderSide(
                        color: Colors.grey,
                      ), // Add border to the button
                    ),
                    child: Image.asset(
                      'assets/google_logo.png',
                      height: 24,
                    ), // Add Google logo in assets
                  ),
                ],
              ),
              SizedBox(height: 15),

              Row(
                mainAxisAlignment: MainAxisAlignment.center,
                children: [
                  TextButton(
                    onPressed: () {},
                    child: Text(
                      "Terms of Service",
                      style: TextStyle(fontSize: 15),
                    ),
                  ),
                  Text(" | "),
                  TextButton(
                    onPressed: () {},
                    child: Text(
                      "Privacy Policy",
                      style: TextStyle(fontSize: 15),
                    ),
                  ),
                ],
              ),
            ],
          ),
        ),
      ),
    );
  }
}
