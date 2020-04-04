/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rappsilber.gui.components.config;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.SwingUtilities;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import rappsilber.gui.SimpleXiGui;
import rappsilber.gui.components.GenericTextPopUpMenu;
import rappsilber.utils.Util;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class TextConfig extends javax.swing.JPanel implements ConfigProvider{

    private BasicConfig basicConfig;
    
    private ArrayList<ActionListener> textConfigListener = new ArrayList<>();
    
    /**
     * Creates new form TextConfig
     */
    public TextConfig() {
        initComponents();
        loadDefaultConfig();
        GenericTextPopUpMenu gtpm = new GenericTextPopUpMenu();
        gtpm.installContextMenu(this);
        
        // enable default load on change
        txtConfig.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                btnDefault.setEnabled(true);
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                btnDefault.setEnabled(true);
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                btnDefault.setEnabled(true);
            }
        });
        
        fbSaveConfig.setLocalPropertyKey("XLink_Config");
        fbSaveConfig.setExtensions(new String[]{".conf",".txt"});
        fbSaveConfig.setDescription("config-files");

    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        fbSaveConfig = new rappsilber.gui.components.FileBrowser();
        btnSave = new javax.swing.JButton();
        btnDefault = new javax.swing.JButton();
        jScrollPane1 = new javax.swing.JScrollPane();
        txtConfig = new javax.swing.JTextArea();
        btnTransfer = new javax.swing.JButton();

        fbSaveConfig.setExtensions(new String[] {"txt"});

        btnSave.setText("Save File");
        btnSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSaveActionPerformed(evt);
            }
        });

        btnDefault.setText("Load Default Config");
        btnDefault.setEnabled(false);
        btnDefault.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnDefaultActionPerformed(evt);
            }
        });

        txtConfig.setColumns(20);
        txtConfig.setRows(5);
        jScrollPane1.setViewportView(txtConfig);

        btnTransfer.setText("To Basic Config");
        btnTransfer.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnTransferActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addComponent(btnDefault)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(layout.createSequentialGroup()
                .addComponent(btnTransfer)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(fbSaveConfig, javax.swing.GroupLayout.DEFAULT_SIZE, 405, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(btnSave))
            .addComponent(jScrollPane1)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(btnDefault)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 214, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(fbSaveConfig, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnSave)
                    .addComponent(btnTransfer))
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    
    public void loadConfig(File f, boolean append) {
        //File f = btnLoadConfig.getFile();
        StringBuffer config = new StringBuffer();
        try {
            if (f!= null) {
                BufferedReader confIn = new BufferedReader(new FileReader(f));
                String line;
                while ((line = confIn.readLine()) != null) {
                    config.append(line);
                    config.append("\n");
                }
                confIn.close();
                if (append)
                    txtConfig.append("\n" + config.toString());
                else 
                    txtConfig.setText(config.toString());
            }
        }catch (IOException ioe) {
            System.err.println(ioe);
        }        
    }
    
    @Override
    public void loadConfig(String config, boolean append, HashSet<String> ignoreSettings) {
        if (ignoreSettings != null)
            for (String i : ignoreSettings) {
                config = config.replaceAll("\n(" + i + ":.*)", "\n# ignored : $1");
            }
        if (append)
            txtConfig.append("\n" + config);
        else 
            txtConfig.setText(config);
    }    
    
    public void appendConfig(String s) {
        txtConfig.append("\n"+s);
    }
    
    public String getConfig() {
        return txtConfig.getText();
    }

    public void setConfig(final String config) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                txtConfig.setText(config);
                txtConfig.setCaretPosition(1);
            }
        });
    }
    
    public void loadDefaultConfig() {
        try {
            BufferedReader br = Util.readFromClassPath(".rappsilber.data.DefaultConfig.conf");
            StringBuffer sb = new StringBuffer();
            String line;
            while ((line = br.readLine()) != null)
                sb.append(line + "\n");
            br.close();
            txtConfig.setText(sb.toString());
            txtConfig.setCaretPosition(1);

        } catch (IOException ex) {
            Logger.getLogger(SimpleXiGui.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    private void btnSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSaveActionPerformed
        FileWriter fw = null;
        try {
            File f = fbSaveConfig.getFile();
            fw = new FileWriter(f);
            fw.write(txtConfig.getText());
            fw.close();
        } catch (IOException ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            }
        }
    }//GEN-LAST:event_btnSaveActionPerformed

    private void btnDefaultActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnDefaultActionPerformed
        loadDefaultConfig();
        btnDefault.setEnabled(false);
    }//GEN-LAST:event_btnDefaultActionPerformed

    private void btnTransferActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnTransferActionPerformed
            
        ActionEvent e = new ActionEvent(this, 0, getConfig());
        for (ActionListener al :textConfigListener) {
            al.actionPerformed(e);
        }
        
    }//GEN-LAST:event_btnTransferActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnDefault;
    private javax.swing.JButton btnSave;
    private javax.swing.JButton btnTransfer;
    private rappsilber.gui.components.FileBrowser fbSaveConfig;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTextArea txtConfig;
    // End of variables declaration//GEN-END:variables

    /**
     * @return the basicConfig
     */
    public BasicConfig getBasicConfig() {
        return basicConfig;
    }

    /**
     * @param basicConfig the basicConfig to set
     */
    public void setBasicConfig(BasicConfig basicConfig) {
        this.basicConfig = basicConfig;
    }

    public void addTransferListener(ActionListener listener) {
        textConfigListener.add(listener);
    }

    public void removeTransferListener(ActionListener listener) {
        textConfigListener.remove(listener);
    }


    
}
