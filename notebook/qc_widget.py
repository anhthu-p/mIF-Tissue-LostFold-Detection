from qtpy.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit
import numpy as np
import zarr
import pandas as pd
import os

class QCWidget(QWidget):

    def __init__(self,
                 viewer,
                 assessor,
                 flagged_df,
                 save_path,
                 resolution_level):
        # what is flagged_df: not_okay biomarker
        # save_path: write df about decision where approve or need to be correct
        super().__init__()
        self.loader            = assessor.loader
        self.biomarker_channel = assessor.biomarker_channel
        self.viewer           = viewer
        self.df = flagged_df.copy().reset_index(drop=True) 
        self.df['decision'] = None
        self.df['notes']    = None
        self.save_path        = save_path        # must be before _load_current
        self.resolution_level = resolution_level     # must be before _load_current
        self.idx              = 0

        self._build_ui()
        self._load_current()


    
    def _build_ui(self):
        # 1. create your 3 labels and store as self.* so other methods can update them
        self.label_name     = QLabel("")
        self.label_flag     = QLabel("")
        self.label_progress = QLabel("")

        #Self note:
        self.notes_input = QLineEdit()
        self.notes_input.setPlaceholderText("Notes (optional)")
        # 2. create approve + reject buttons with shortcuts A and R
        button_yes = QPushButton("Approve [A]")
        button_no = QPushButton("Reject [R]")
        button_yes.setShortcut('A')
        button_no.setShortcut('R')

        # If go back:
        button_back = QPushButton("Back [B]")
        button_back.setShortcut('B')
        button_back.clicked.connect(self._go_back)

        # 3. connect buttons to self._record
        button_yes.clicked.connect(self.on_approve)
        button_no.clicked.connect(self.on_reject)
        # 4. put buttons side by side in QHBoxLayout
        btn_row = QHBoxLayout()
        btn_row.addWidget(button_yes)
        btn_row.addWidget(button_no)
        btn_row.addWidget(button_back)
        
        # 5. stack everything in QVBoxLayout
        layout = QVBoxLayout()
        layout.addWidget(self.label_name)
        layout.addWidget(self.label_flag)
        layout.addWidget(self.label_progress)
        layout.addWidget(self.notes_input)
        layout.addLayout(btn_row)   # note: addLayout not addWidget for nested layouts
        # 6. self.setLayout(...)
        self.setLayout(layout)
    def _load_current(self):
        if self.idx >= len(self.df):
            self.label_name.setText("Done!")
            self.label_flag.setText("")
            self.label_progress.setText(f"{len(self.df)} / {len(self.df)}")
            return
        self._update_labels()
        self._load_image()

    def _update_labels(self):
        row = self.df.iloc[self.idx]
        self.label_name.setText(row['biomarker'])
        self.label_flag.setText(row['flag'])
        self.label_progress.setText(f"{self.idx + 1} / {len(self.df)}")
        existing_note = self.df.loc[self.idx, 'notes']
        self.notes_input.setText(str(existing_note) if existing_note else "")

    def _load_image(self):
        row      = self.df.iloc[self.idx]
        img_full = zarr.open(self.loader.image_store, mode='r')
        self.viewer.layers.clear()

        img = img_full[str(self.resolution_level)][int(row['channel_index'])]
        low, high = np.percentile(img, (2.5, 99.8))
        self.viewer.add_image(img, name=row['biomarker'],
                            contrast_limits=[low, high],
                            colormap='gray', blending='additive')

        dna_idx = int(self.biomarker_channel['dna'])
        img_dna = img_full[str(self.resolution_level)][dna_idx]
        low, high = np.percentile(img_dna, (2.5, 99.8))
        self.viewer.add_image(img_dna, name='dna',
                            contrast_limits=[low, high],
                            colormap='blue', blending='additive')

    
    def on_approve(self):
        self._record('approved')
    def on_reject(self):
        self._record('rejected')
    def _go_back(self):
        if self.idx > 0:
            self.idx -= 1
            self.df.loc[self.idx, 'decision'] = None  # clear previous decision
            self.df.loc[self.idx, 'notes']    = None
            self._load_current()


    def _record(self, decision):
        self.df.loc[self.idx, 'decision'] = decision
        self.df.loc[self.idx, 'notes']    = self.notes_input.text()
        self.notes_input.clear()
        save_dir = os.path.dirname(self.save_path)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
        self.df.to_csv(self.save_path, index=False)
        self.idx += 1
        self._load_current()
