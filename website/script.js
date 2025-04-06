// Wait for the DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function () {
    // ===== SMOOTH SCROLLING NAVIGATION =====
    const navLinks = document.querySelectorAll('nav a');

    navLinks.forEach(link => {
        link.addEventListener('click', function (e) {
            const href = this.getAttribute('href');

            // Only prevent default behavior for internal links (those starting with '#')
            if (href.startsWith('#')) {
                e.preventDefault();

                // Remove active class from all links
                navLinks.forEach(link => link.classList.remove('active'));

                // Add active class to clicked link
                this.classList.add('active');

                // Get the target section id from href attribute
                const targetId = href.substring(1);
                const targetSection = document.getElementById(targetId);

                if (targetSection) {
                    // Smooth scroll to target section
                    targetSection.scrollIntoView({
                        behavior: 'smooth'
                    });
                }
            }
        });
    });

    // ===== LEARN MORE BUTTON SMOOTH SCROLL =====
    const learnMoreButton = document.querySelector('.hero-buttons .button');
    if (learnMoreButton) {
        learnMoreButton.addEventListener('click', function (e) {
            e.preventDefault();
            const targetId = this.getAttribute('href').substring(1);
            const targetSection = document.getElementById(targetId);
            if (targetSection) {
                targetSection.scrollIntoView({
                    behavior: 'smooth'
                });
            }
        });
    }

    // ===== CONTACT FORM FUNCTIONALITY =====
    const contactForm = document.querySelector('.contact-form form');
    
    if (contactForm) {
        // Form elements
        const nameInput = document.getElementById('name');
        const emailInput = document.getElementById('email');
        const messageInput = document.getElementById('message');
        const submitButton = document.querySelector('.contact-form button');

        // Create notification element
        const notification = document.createElement('div');
        notification.className = 'notification';
        notification.style.cssText = `
            display: none;
            background-color: #4CAF50;
            color: white;
            padding: 15px;
            margin-top: 20px;
            border-radius: 8px;
            text-align: center;
            transition: all 0.5s ease;
            opacity: 0;
            transform: translateY(-20px);
        `;
        contactForm.appendChild(notification);

        // Form validation functions
        function validateForm() {
            return nameInput.value.trim() !== '' &&
                emailInput.value.trim() !== '' &&
                messageInput.value.trim() !== '';
        }

        function validateEmail(email) {
            const re = /^(([^<>()\[\]\\.,;:\s@"]+(\.[^<>()\[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;
            return re.test(String(email).toLowerCase());
        }

        // Form field listeners
        [nameInput, emailInput, messageInput].forEach(input => {
            input.addEventListener('input', checkFormValidity);
        });

        // Update button state based on form validity
        function checkFormValidity() {
            const isValid = validateForm();
            submitButton.disabled = !isValid;
            submitButton.style.opacity = isValid ? '1' : '0.6';
            submitButton.style.cursor = isValid ? 'pointer' : 'not-allowed';
        }

        // Initial button state
        submitButton.disabled = true;
        submitButton.style.opacity = '0.6';
        submitButton.style.cursor = 'not-allowed';

        // Form submission handler
        contactForm.addEventListener('submit', async function (event) {
            event.preventDefault();
            
            if (!validateForm()) return;
            if (!validateEmail(emailInput.value.trim())) {
                showNotification('Please enter a valid email address.', 'error');
                return;
            }

            // Update button state
            submitButton.textContent = 'Sending...';
            submitButton.disabled = true;
            submitButton.style.opacity = '0.6';

            // Prepare form data
            const formData = {
                name: nameInput.value.trim(),
                email: emailInput.value.trim(),
                message: messageInput.value.trim()
            };

            try {
                const response = await sendEmail(formData);
                submitButton.textContent = 'Send Message';

                if (response.success) {
                    showNotification('Message sent successfully!', 'success');
                    contactForm.reset();
                    checkFormValidity();
                } else {
                    throw new Error(response.message);
                }
            } catch (error) {
                submitButton.textContent = 'Send Message';
                submitButton.disabled = false;
                submitButton.style.opacity = '1';
                submitButton.style.cursor = 'pointer';
                showNotification('Failed to send message: ' + error.message, 'error');
            }
        });

        // Notification display function
        function showNotification(message, type) {
            notification.textContent = message;
            notification.style.backgroundColor = type === 'success' ? '#4CAF50' : '#f44336';
            notification.style.display = 'block';

            setTimeout(() => {
                notification.style.opacity = '1';
                notification.style.transform = 'translateY(0)';
            }, 10);

            setTimeout(() => {
                notification.style.opacity = '0';
                notification.style.transform = 'translateY(-20px)';
                setTimeout(() => {
                    notification.style.display = 'none';
                }, 500);
            }, 3000);
        }
    }

    // ===== SCROLL-BASED NAVIGATION HIGHLIGHTING =====
    window.addEventListener('scroll', function () {
        const scrollPosition = window.scrollY;
        const sections = document.querySelectorAll('section[id]');

        sections.forEach(section => {
            const sectionTop = section.offsetTop - 100;
            const sectionHeight = section.offsetHeight;
            const sectionId = section.getAttribute('id');

            if (scrollPosition >= sectionTop && scrollPosition < sectionTop + sectionHeight) {
                navLinks.forEach(link => link.classList.remove('active'));
                const activeLink = document.querySelector(`nav a[href="#${sectionId}"]`);
                if (activeLink) {
                    activeLink.classList.add('active');
                }
            }
        });
    });
});
// No sendEmail function needed - Netlify handles it automatically
document.querySelector('form[name="contact"]').addEventListener('submit', async function(e) {
    e.preventDefault();
    
    try {
        const formData = new FormData(this);
        await fetch('/', {
            method: 'POST',
            headers: { "Content-Type": "application/x-www-form-urlencoded" },
            body: new URLSearchParams(formData).toString()
        });
        
        showNotification('Message sent successfully!', 'success');
        this.reset();
    } catch (error) {
        showNotification('Failed to send message: ' + error.message, 'error');
    }
});
