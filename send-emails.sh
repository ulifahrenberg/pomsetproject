#!/bin/bash

#ADDRS=private/emails-test.txt
ADDRS=private/emails.txt

FROM=fahrenberg@gmail.com

EMAILCONT=announce.txt

#SUBJECT='(i)Po(m)set Project Online Seminar: Safa Zouari tomorrow 17 March 11:00 CET'
#SUBJECT='(i)Po(m)set Project Online Seminar: Sergio Rajsbaum Friday 14 April'
#SUBJECT='(i)Po(m)set Project Online Seminar: Mikołaj Bojańczyk tomorrow 15 December 10:30 CET 🎄'
SUBJECT='Third International Workshop on Pomsets and Related Structures'

#TMPDIR=$(mktemp -d)
TMPDIR=/tmp

echo "Sending emails to all people in $ADDRS"

for LINE in $(cat $ADDRS); do
    #echo $LINE
    FIRST=$(echo $LINE|awk -F, '{print $2}')
    EMAIL=$(echo $LINE|awk -F, '{print $3}')
    echo "Sending to $FIRST at $EMAIL"
    echo "From: $FROM" > $TMPDIR/email.txt
    echo "To: $EMAIL" >> $TMPDIR/email.txt
    echo "Subject: $SUBJECT" >> $TMPDIR/email.txt
    echo "" >> $TMPDIR/email.txt
    #echo "Dear $FIRST," >> $TMPDIR/email.txt
    echo "Dear colleague," >> $TMPDIR/email.txt
    echo "" >> $TMPDIR/email.txt
    cat $EMAILCONT >> $TMPDIR/email.txt
    cat $TMPDIR/email.txt | ssmtp $EMAIL
done

echo 'Done'
