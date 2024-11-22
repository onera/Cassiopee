# Send a notification by email to a list of recipients
# Please set the environment variable CASSIOPEE_EMAIL
# Usage: python notify.py --recipients='a.b@onera.fr c.d@onera.fr' --subject='subject of the message' --message='a message to send'
import os

# Send an email
def notify(sender=None, recipients=[], messageSubject="", messageText=""):
    try:
        import smtplib
        from email.mime.text import MIMEText
        from email.mime.multipart import MIMEMultipart

        if sender is None:
            if os.getenv('CASSIOPEE_EMAIL') is None:
                if os.getenv('USER') is None:
                    print("Sender email address not found.")
                    return
                else: sender = os.getenv('USER')+'@onera.fr'
            else: sender = os.getenv('CASSIOPEE_EMAIL')
        if isinstance(recipients, str): recipients = [recipients]
        if not recipients: recipients = [sender]

        msg = MIMEMultipart()
        msg['Subject'] = messageSubject
        msg['From'] = sender
        msg['To'] = ", ".join(recipients)
        msg['Cc'] = sender
        msg.preamble = 'Sent by Cassiopee.'
        if messageText:
            msg.attach(MIMEText(messageText, 'plain'))
        s = smtplib.SMTP('localhost')
        s.sendmail(sender, recipients, msg.as_string())
        s.quit()
    except: return

# Parse command-line arguments
def parseArgs():
    import argparse
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--message", type=str, default='',
                        help="Single-quoted message")
    parser.add_argument("-r", "--recipients", type=str, default='',
                        help="Single-quoted space-separated list of recipients")
    parser.add_argument("-s", "--subject", type=str, default='',
                        help="Single-quoted message subject")
    # Parse arguments
    return parser.parse_args()

# Main
if __name__ == '__main__':
    script_args = parseArgs()
    recipients = script_args.recipients.split(' ')
    if not recipients[0]: recipients = []
    notify(recipients=recipients, messageSubject=script_args.subject,
           messageText=script_args.message)
