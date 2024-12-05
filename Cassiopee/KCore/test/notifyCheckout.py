# Notify a user about the checkout status of different prods, either in the
# terminal window (default) or by email (set the environment variable
# CASSIOPEE_EMAIL and run on localhost)
# Usage: python notifyCheckout.py
#        python notifyCheckout.py --email --recipients='a.b@onera.fr c.d@onera.fr'
import sys
from time import strptime, strftime

# Parse command-line arguments
def parseArgs():
    import argparse
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--email", action="store_true",
                        help="Email results. Default: print in terminal")
    parser.add_argument("-r", "--recipients", type=str, default='',
                        help="Single-quoted space-separated list of recipients")
    # Parse arguments
    return parser.parse_args()

# Main
if __name__ == '__main__':
    try:
        import KCore.Dist as Dist
        from KCore.notify import notify
    except ImportError:
        print("Error: KCore is required to execute notifyCheckout.py")
        sys.exit()

    script_args = parseArgs()
    recipients = script_args.recipients.split(' ')
    if not recipients[0]: recipients = []

    # Check checkout status
    log_entries = []
    with open('/stck/cassiope/git/logs/checkout_status.txt', 'r') as f:
        for line in f:
            log_entries.append(line.strip().split(' - '))
    log_entries.sort(key=lambda x: x[1], reverse=True)

    # Do not send a notification when everything is OK
    if not any('FAILED' in log_machine for log_machine in log_entries):
        if not script_args.email: print("[Checkout Cassiopee] State: OK")
        sys.exit()

    # Get git info
    cassiopeeIncDir = '/stck/cassiope/git/Cassiopee/Cassiopee'
    gitOrigin = Dist.getGitOrigin(cassiopeeIncDir)
    gitBranch = Dist.getGitBranch(cassiopeeIncDir)
    gitHash = Dist.getGitHash(cassiopeeIncDir)[:7]
    gitInfo = "Git origin: {}\nGit branch: {}\nCommit hash: {}".format(
        gitOrigin, gitBranch, gitHash)

    messageSubject = "[Checkout Cassiopee] State: FAILED"
    messageText = "Pulling updates for Cassiopee, Fast and all "\
        "PModules:\n{}\n\n{}\n\n".format(52*'-', gitInfo)
    messageText += '{:^20} | {:^15} | {:^30} | {:^10}\n{}\n'.format(
        "PROD.", "PCKG.", "DATE", "STATUS", 85*'-')
    for log_machine in log_entries:
        prod = log_machine[0]
        pckg = log_machine[1]
        date = strptime(log_machine[2], "%y%m%d-%H%M%S")
        date = strftime("%d/%m/%y at %T", date)
        status = log_machine[3]
        messageText += '{:^20} | {:^15} | {:^30} | {:^10}\n'.format(
            prod, pckg, date, status)

    messageText += '\n\nIf the prod. you wish to use is marked as FAILED, '\
        'please contact the maintainers:\nchristophe.benoit@onera.fr, '\
        'vincent.casseau@onera.fr'

    if script_args.email:
        notify(recipients=recipients,
               messageSubject=messageSubject,
               messageText=messageText)
    else:
        print("{0}\n|{1:^65}|\n{0}\n{2}".format(67*'-', messageSubject, messageText))
