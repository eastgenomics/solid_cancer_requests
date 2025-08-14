"""
Scraper to pull driver mutation data from cbioportal.

Requires setting a location to download files to (`download_dir`) and path
to chrome driver (`chrome_driver`) downloaded from https://googlechromelabs.github.io/chrome-for-testing/.
The chrome driver version needs to match the version of chrome installed.

Genes to download data for should be provided as input in a single column file
with one gene symbol per row.
"""

from pathlib import Path
import sys
import time

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import undetected_chromedriver as uc

download_dir = "/home/jethro/Projects/solid_cancer_requests/files/"
chrome_driver = "/home/jethro/Projects/solid_cancer_requests/chromedriver-linux64/chromedriver"

with open(sys.argv[1]) as fh:
    genes = fh.read().splitlines()

options = uc.ChromeOptions()
options.add_experimental_option(
    "prefs",
    {"download.default_directory": download_dir},
)
options.add_argument("--no-sandbox")
options.add_argument("--disable-dev-shm-usage")

service = Service(chrome_driver)
driver = webdriver.Chrome(service=service, options=options)
driver = uc.Chrome(options=options)

driver.get("https://www.cbioportal.org/")

print("Loading page for login")
login_url = "https://genie.cbioportal.org"
driver.get(login_url)

# need to allow time to first manually log into Google account,
# only needs doing once
input("Press any key after logging in and page has loaded")

failed = []

for idx, gene in enumerate(genes, 1):
    print(f"[{idx}/{len(genes)}] Loading page for {gene}")

    try:
        target_url = (
            "https://genie.cbioportal.org/results/mutations?cancer_study_list=genie_public&"
            "Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants"
            f"%2Ccna&case_set_id=genie_public_cnaseq&gene_list={gene}&geneset_list=%20&tab_index=tab_visualize&Action=Submit"
        )

        driver.get(target_url)

        wait = WebDriverWait(
            driver, 300
        )  # set long timeout since page is slow to load

        # find the download button in the page
        download_button = wait.until(
            EC.element_to_be_clickable(
                (
                    By.XPATH,
                    "//button[@aria-label='Download TSV']",
                )
            )
        )
        download_button.click()

        print(f"Downloading {gene} data, waiting until complete...")

        downloaded_file = Path(Path.cwd() / "files/table.tsv")

        while not downloaded_file.exists():
            time.sleep(1)

        rename_file = Path(downloaded_file.parent, f"{gene}.tsv")
        downloaded_file.rename(rename_file)

        print(f"Downloaded data to {str(rename_file)}")
    except Exception as err:
        print(f"Error getting data for {gene}: {err}")
        failed.append(f"{gene} - {err}")

print(f"Downloaded data for {len(genes)} genes")

if failed:
    failed = "\n\t".join(failed)
    print(f"Warning: failed to get data for {len(genes)} genes.\n{failed}")

driver.quit()
