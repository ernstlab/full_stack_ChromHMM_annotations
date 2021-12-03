# downloaded from https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655447/Programming+Interface on May 01 2021
#!/usr/bin/env python

import sys;
from threading import Thread;
import urllib.request;
import time;

if len(sys.argv) != 2:
	sys.stderr.write("Usage: %s <joblist>\n" % (sys.argv[0]));
	sys.exit(1);

NUM_REQUESTS = 5 # hard limit by great server, raising this will not get more jobs since the server will deny requests
JOBLIST = open(sys.argv[1]).readlines();

def main():
	for r in range(NUM_REQUESTS):
		t = Agent(JOBLIST[r::NUM_REQUESTS]); # split job list into the number requests possible
		t.start();

class Agent(Thread):
	def __init__(self, joblist):
		Thread.__init__(self);
		self.jobs = [];
		for job in joblist:
			self.jobs.extend(list(map(lambda x: x.strip().split(), joblist)));
	
	def run(self):
		jobs = self.jobs;
		for output, url in jobs:
			req = urllib.request.Request(url);
			retry = True;
			while retry:
				try:
					retry = False;
					result = urllib.request.urlopen(req);
				except urllib.request.HTTPError as error: 
					if(error.getcode() != 500):
						sys.stderr.write("Expected Error: [HTTP Error 500: Internal Server Error]. Actual Error: [%s].\n" % (error));
						sys.stderr.write("\tQuitting due to unexpected error.\n");
						sys.exit(1);
					#sys.stderr.write("[%s] for output [%s]. Server likely busy. Will retry.\n" % (error, output));	# Remove comment if you want to see retries
					retry = True;
					time.sleep(10); # if request gets denied, wait 10 seconds and try again
				if(not retry):
					f = open(output, 'wb');
					f.write(result.read());
					f.close();
					print ("STATUS: Done, REQUEST: %s, OUTPUT: %s" % (url, output));
if __name__ == '__main__':
	main();
