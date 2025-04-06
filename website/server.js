const express = require('express');
const nodemailer = require('nodemailer');
const cors = require('cors');
const path = require('path');
require('dotenv').config();

const app = express();

// Configure CORS more explicitly
app.use(cors({
    origin: ['http://localhost:5500', 'http://127.0.0.1:5500'], // Allow both localhost variations
    methods: ['GET', 'POST'],
    credentials: true
}));

// Middleware to parse JSON bodies
app.use(express.json());

// Serve static files (including script.js)
app.use(express.static(path.join(__dirname)));

// Email configuration
const EMAIL_HOST = 'smtp.gmail.com';
const EMAIL_PORT = 587;
const EMAIL_HOST_USER = process.env.EMAIL_HOST_USER || 'medixir.ai@gmail.com';
const EMAIL_HOST_PASSWORD = process.env.EMAIL_HOST_PASSWORD || 'ngvz loba obfg kdav';
const RECIPIENT_EMAIL = process.env.RECIPIENT_EMAIL || 'medixir.ai@gmail.com';

// Create email transporter
const transporter = nodemailer.createTransport({
    host: EMAIL_HOST,
    port: EMAIL_PORT,
    secure: false,
    auth: {
        user: EMAIL_HOST_USER,
        pass: EMAIL_HOST_PASSWORD
    }
});

// Route to serve the main index.html file
app.get('/', (req, res) => {
    res.sendFile(path.join(__dirname, 'index.html'));
});

// Route to send email
app.post('/send-email', async (req, res) => {
    try {
        const { name, email, message } = req.body;
        
        // Log the received data for debugging
        console.log(`Received form data: Name=${name}, Email=${email}, Message=${message}`);
        
        // Create email options
        const mailOptions = {
            from: EMAIL_HOST_USER,
            to: RECIPIENT_EMAIL,
            subject: `New message from ${name} via AI DrugLab`,
            text: `
            You have received a new message from your website:
            
            Name: ${name}
            Email: ${email}
            
            Message:
            ${message}
            `
        };
        
        console.log("Attempting to send email...");
        
        try {
            await transporter.sendMail(mailOptions);
            console.log("Email sent successfully!");
            res.json({ success: true, message: "Email sent successfully" });
        } catch (error) {
            console.error("SMTP Error:", error);
            res.status(500).json({ 
                success: false, 
                message: `SMTP Error: ${error.message}`
            });
        }
        
    } catch (error) {
        console.error("Error sending email:", error);
        res.status(500).json({ 
            success: false, 
            message: error.message 
        });
    }
});

// Test route to check email configuration
app.get('/test-email', async (req, res) => {
    try {
        await transporter.verify();
        res.json({ 
            success: true, 
            message: "SMTP connection test successful" 
        });
    } catch (error) {
        res.status(500).json({ 
            success: false, 
            message: `SMTP connection test failed: ${error.message}` 
        });
    }
});

const PORT = 5500;
const HOST = '127.0.0.1';

// Start the server
if (EMAIL_HOST_PASSWORD === "ngvz loba obfg kdav") {
    console.log("=".repeat(50));
    console.log("WARNING: You need to replace 'ngvz loba obfg kdav' with your actual App Password");
    console.log("Follow the instructions in the comments to set up an App Password");
    console.log("=".repeat(50));
}

app.listen(PORT, HOST, () => {
    console.log("=".repeat(50));
    console.log("Server is running!");
    console.log(`Access your website at: http://${HOST}:${PORT}/`);
    console.log(`Visit http://${HOST}:${PORT}/test-email to test email connectivity`);
    console.log("=".repeat(50));
});