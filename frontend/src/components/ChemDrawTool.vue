<template>
  <div :id="jsmeContainerName"></div>
  <input :id="smilesContainerName" hidden="hidden" v-model="smiles">
</template>

<script>
import { defineComponent } from 'vue'
export default defineComponent(
  {
    name: 'ChemDrawTool',
    emits: ['smiles'],
    data () {
      return {
        smiles: '',
        jsmeContainerName: 'jsmeContainer',
        smilesContainerName: 'smilesContainer',
        storeSmilesFunctionName: 'storeSmiles',
        jsmeScriptName: 'jsmeScript',
        jsmeOnLoadScriptName: 'jsmeOnLoad'
      }
    },
    mounted () {
      const jsmeScript = document.createElement('script')
      jsmeScript.src = 'jsme/jsme.nocache.js'
      jsmeScript.type = 'text/javascript'
      jsmeScript.language = 'javascript'
      jsmeScript.id = this.jsmeScriptName
      document.head.appendChild(jsmeScript)

      const jsmeOnLoad = document.createElement('script')
      jsmeOnLoad.id = this.jsmeOnLoadScriptName
      jsmeOnLoad.innerHTML = `
        function jsmeOnLoad() {
            jsmeApplet = new JSApplet.JSME("${this.jsmeContainerName}", "380px", "340px");
            jsmeApplet.setCallBack("AfterStructureModified", ${this.storeSmilesFunctionName});
        }
      `
      document.head.appendChild(jsmeOnLoad)

      const storeSmiles = document.createElement('script')
      storeSmiles.id = this.storeSmilesFunctionName
      storeSmiles.innerHTML = `
        function ${this.storeSmilesFunctionName}(event) {
          smiles = event.src.smiles();
          const smilesContainer = document.getElementById("${this.smilesContainerName}");
          smilesContainer.value = smiles;
          smilesContainer.dispatchEvent(new Event('input'));
        }
      `
      document.head.appendChild(storeSmiles)
    },
    unmounted () {
      // Clean up jsme editor on route change
      document.getElementById(this.jsmeScriptName).remove()
      document.getElementById(this.jsmeOnLoadScriptName).remove()
      document.getElementById(this.storeSmilesFunctionName).remove()
    },
    watch: {
      smiles (value) {
        this.smiles = value
        this.$emit('smiles', this.smiles)
      }
    }
  }
)
</script>
